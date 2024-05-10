import bpy
import itertools
from mathutils import Matrix, Quaternion, Vector
import multiprocessing as mp
import os
from pathlib import Path


z_to_y = Matrix([
    [ 1, 0, 0 ],
    [ 0, 0, 1 ],
    [ 0, -1, 0 ]]).to_quaternion()

def transpose(lists):
    transposed = [[] for _ in lists[0]]
    for list in lists:
        for i, l in enumerate(list):
            transposed[i].append(l)
    return transposed

def quaternion_product(q1, q2):
    q3 = [0, 0, 0, 0]
    q3[0] = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z
    q3[1] = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y
    q3[2] = q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x
    q3[3] = q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w
    return Quaternion(q3)


class Vertex:
    def __init__(self, index):
        self.index = index
        self.tail_position = Vector((0, 0, 0))
        self.radius = 0
        self.on_leaf = False
        self.fixed = False
        
    def __str__(self):
        x = self.tail_position.x
        y = self.tail_position.y
        z = self.tail_position.z
        on_leaf = 1 if self.on_leaf else 0
        fixed = 1 if self.fixed else 0
        return f'{self.index},{x},{y},{z},{self.radius},{on_leaf},{fixed}'



def get_bone_radius(bone):
    return (bone.head_radius + bone.tail_radius) / 2

def is_leaf_like(bone, depth):
    if len(bone.children) == 0:
        return True
    if depth == 1:
        return False
    for child in bone.children:
        if is_leaf_like(child, depth - 1):
            return True
        else:
            return False

class PlantSimulation:

# Exporting the plant as an armature
    def __init__(self, obj):
        self.global_scale = 10
        assert(obj.type == 'ARMATURE')
        self.obj = obj
        self.armature = obj.data
        self.pose = obj.pose
        self.bone_count = len(self.armature.bones)
        self.leaf_depth = 3

        # Make a reverse index for bones
        self.reverse_index = {}
        for i, bone in enumerate(self.armature.bones):
            self.reverse_index[bone] = i
        for i, bone in enumerate(self.pose.bones):
            self.reverse_index[bone] = i
        
        # Compute rest index
        # Converts y-up coordinate to world coordinate
        self.rest_quaternions = [b.matrix_local.to_quaternion() for b in self.armature.bones]


    def export(self, out_file_path):
        
        # Make a list of segments
        self.vertices = []
        for i,bone in enumerate(self.armature.bones):
            assert(get_bone_radius(bone) > 0)
            assert(bone.head != bone.tail)
            new_vertex = Vertex(i)
            new_vertex.tail_position = bone.tail_local / self.global_scale
            new_vertex.radius = get_bone_radius(bone) / self.global_scale
            # Arbitrary.
            new_vertex.on_leaf = is_leaf_like(bone, self.leaf_depth)
            new_vertex.fixed = bone.parent is None
            self.vertices.append(new_vertex)

        # Make a list of edges
        self.edges = []
        for i, bone in enumerate(self.armature.bones):
            if bone.parent is None:
                continue
            parent_index = self.reverse_index[bone.parent]
            self.edges.append((parent_index, i))

        # Export a plant file
        with open(out_file_path, 'w') as out_file:
            print(f'verts {len(self.vertices)}', file=out_file)
            for vertex in self.vertices:
                print(str(vertex), file=out_file)
                
            print(f'edges {len(self.edges)}', file=out_file)
            for edge in self.edges:
                print(f'{edge[0]},{edge[1]}', file=out_file)
         
         
    def read_pose_file(self, in_file_path):

        with open(in_file_path) as in_file:
            line = in_file.readline().split(' ')
            assert(line[0] == 'rods')
            assert(int(line[1]) == self.bone_count)

            # Read the length and rotations of each bone
            # Remember to use self.global_scale
            poses = [PlantSimulation.read_line(l) for l in in_file]
            
        return poses
    
    def read_line(line):
        tokens = line.split(',')
        length = float(tokens[1])
        w = float(tokens[2])
        x = float(tokens[3])
        y = float(tokens[4])
        z = float(tokens[5])
        q = Quaternion((w, x, y, z))
        return q, length
    
    def compute_initial_poses(self, in_file_path):
        self.initial_poses = self.read_pose_file(in_file_path)
        
    def read_pose(self, path):
        poses = self.read_pose_file(path)
        
        for i, bone in enumerate(self.pose.bones):
            q, length = poses[i]
            q_rest = self.rest_quaternions[i]
            q_initial, l_initial = self.initial_poses[i]
            if bone.parent is None:
                q_parent = Quaternion((1, 0, 0, 0))
                q_parent_initial = Quaternion((1, 0, 0, 0))
                q_parent_rest = z_to_y.conjugated()
            else:
                q_parent, _ = poses[self.reverse_index[bone.parent]]
                q_parent_initial, _ = self.initial_poses[self.reverse_index[bone.parent]]
                q_parent_rest = self.rest_quaternions[self.reverse_index[bone.parent]]
            bone.rotation_mode = 'QUATERNION'


            q_from_parent = q_parent.conjugated() @ q
            q_from_parent_initial = q_parent_initial.conjugated() @ q_initial
            q_from_parent_rest = q_parent_rest.conjugated() @ q_rest


            # There are actually four spaces:
            # Blender global space (automatic from bone local space)
            # Bone parent space
            # Blender bone local space
            # Physics local space (relative to parent, Easily calculated from physics local space)
            # Physics bone global space (from physics simulation)
            
            #arbitrary = Quaternion((1, 2, 3, 4)).normalized()
            #q = arbitrary @ Matrix.Rotation(0.1, 4, 'X').to_quaternion()
            #q_parent = arbitrary @ Quaternion((1, 0, 0, 0))
            

            physics_parent_to_global = q_parent
            physics_global_to_parent = physics_parent_to_global.conjugated() # Trivial
            q_physics_from_parent = physics_global_to_parent @ q
            # If we apply q_physics_from_parent to q_physics_local, we should get q_parent
            # the parent-to-local transformations are invisible transformations that are fixed in edit mode
            # Remember that q_initial is global. (Shouldn't actually matter for linear transformations.)
            physics_local_to_parent = q_parent_initial.conjugated() @ q_initial
            physics_parent_to_local = physics_local_to_parent.conjugated() # Trivial
            q_physics_local = physics_parent_to_local @ q_physics_from_parent
            # physics_local_to_parent @ local_local_q should equal q_from_parent
            
            # Experimentation required
            physics_to_blender = Matrix([
                [ 1, 0, 0 ],
                [ 0, 0, 1 ],
                [ 0, -1, 0 ],
            ]).to_quaternion()

            q_blender_from_parent = physics_to_blender @ q_physics_from_parent @ physics_to_blender.conjugated()
            
            # Fixed for any bone during edit mode
            blender_local_to_parent = q_parent_rest.conjugated() @ q_rest
            blender_parent_to_local = blender_local_to_parent.conjugated() # Trivial
            
            q_blender_local = blender_parent_to_local @ q_blender_from_parent
            
            bone.rotation_quaternion = q_blender_local
            
            bone.scale = Vector((1.0, 1.0, length / l_initial))
        
    def import_animations(self, physics_path, max_frames=-1):
        self.compute_initial_poses(physics_path / 'frame-00000.txt')
                
        frame_files = os.listdir(physics_path)
        frame_files.sort()
        if max_frames != -1:
            frame_files = frame_files[max_frames]
        frame_paths = [physics_path / file for file in frame_files]
        
        # Delete all previous keyframes
        for fcurve in self.obj.animation_data.action.fcurves:
            fcurve.keyframe_points.clear()
        self.obj.animation_data.action.fcurves.clear()
        
        print('Reading pose files')
        
        for frame, path in enumerate(frame_paths):
            self.read_pose(path)
            
            for i, bone in enumerate(self.pose.bones):
                bone.keyframe_insert(data_path='rotation_quaternion', frame=frame)
                bone.keyframe_insert(data_path='scale', frame=frame)

        print('Done importing animations')

      
# Select the armature first!

obj = bpy.context.object
assert(obj.type == 'ARMATURE')

simulation = PlantSimulation(obj)

simulation.export('exported-plant.txt')

# After this step, run run-simulation to generate animation data.


base_path = Path('frames/blender')
physics_path = base_path / 'physics'
water_path = base_path / 'water'


simulation.compute_initial_poses(physics_path / 'frame-00000.txt')
#simulation.read_pose(physics_path / 'frame-00000.txt')
simulation.import_animations(physics_path)