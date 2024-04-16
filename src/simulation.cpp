#include "simulation.h"
#include "graphics/meshloader.h"

#include <iostream>

using namespace Eigen;

#define NUM_VERTEX 16
#define NUM_ITER 1
#define C_PER_VERTEX 8

Simulation::Simulation() {}

std::vector<Vector3d> vertices;

// each vertex has three degree of freedom
Eigen::MatrixXd M_inverse = Eigen::MatrixXd::Ones(3*NUM_VERTEX, 3*NUM_VERTEX);

const double alpha = 0.5; //stiffness, assume the same for all vertcies now

Eigen::MatrixXd x = Eigen::MatrixXd::Zero(3, NUM_VERTEX);
Eigen::MatrixXd v = Eigen::MatrixXd::Zero(3, NUM_VERTEX);
Eigen::MatrixXd a = Eigen::MatrixXd::Zero(3, NUM_VERTEX);

Eigen::MatrixXd x0 = Eigen::MatrixXd::Zero(3, NUM_VERTEX);

int findRelateVertex(int currVertex, int c);
double evalConstraint(int currVertex, int c);
double evalScalar(int currVertex, int j, int c);
Eigen::VectorXd evalGradient(int currVertex, int j, int c);

void Simulation::init()
{
    //6  4  2  0
    //7  5  3  1
    //11 10 9  8
    //15 14 13 12
    std::vector<Vector4i> tets;
    std::vector<Vector3i> faces;
    vertices.emplace_back(0,0,0);
    vertices.emplace_back(0,-1,0);

    vertices.emplace_back(1,0,0);
    vertices.emplace_back(1,-1,0);
    faces.emplace_back(1, 0, 2);
    faces.emplace_back(1, 2, 3);
    vertices.emplace_back(2,0,0);
    vertices.emplace_back(2,-1,0);
    faces.emplace_back(3, 2, 4);
    faces.emplace_back(3, 4, 5);
    vertices.emplace_back(3,0,0);
    vertices.emplace_back(3,-1,0);
    faces.emplace_back(5, 4, 6);
    faces.emplace_back(5, 6, 7);

    vertices.emplace_back(0,-2,0);
    vertices.emplace_back(1,-2,0);
    faces.emplace_back(8, 1, 2);
    faces.emplace_back(8, 2, 9);
    vertices.emplace_back(2,-2,0);
    faces.emplace_back(9, 3, 5);
    faces.emplace_back(9, 5, 10);
    vertices.emplace_back(3,-2,0);
    faces.emplace_back(10, 5, 7);
    faces.emplace_back(10, 7, 11);

    vertices.emplace_back(0,-3,0);
    vertices.emplace_back(1,-3,0);
    faces.emplace_back(12, 8, 9);
    faces.emplace_back(12, 9, 13);
    vertices.emplace_back(2,-3,0);
    faces.emplace_back(13, 9, 10);
    faces.emplace_back(13, 10, 14);
    vertices.emplace_back(3,-3,0);
    faces.emplace_back(14, 10, 11);
    faces.emplace_back(14, 11, 15);

    for(int i = 0; i < NUM_VERTEX; i++){
        x.col(i) = vertices[i];
        v.col(i) = Vector3d(0.0, 0.0, 0.0);
        if(i != 0 || i != 6){
            a.col(i) = Vector3d(0.0, -1.0, 0.0);
        }
        x0.col(i) = vertices[i];
    }

    m_shape.init(vertices, faces, tets);
    m_shape.setModelMatrix(Affine3f(Eigen::Translation3f(0, 2, 10)));

    // initGround();
}

void Simulation::update(double seconds)
{
    Eigen::MatrixXd x_tilde = x + seconds * v + seconds * seconds * a;
    Eigen::MatrixXd lambda = Eigen::MatrixXd::Zero(1, NUM_VERTEX * C_PER_VERTEX);
    // TODO: replace itr < NUMITR with  L-Infinity Norm
    for(int itr = 0 ; itr < NUM_ITER; itr++){

        //for all constraint， iteratively solve each constraint -> for large plant, replace this with a direct solver
        //  for the cloth, each vertex has eight constaint: horizontal, vertical and diagonal neighbors:
        //  example constraint C(v5, v4) = dist(v5, v4) - initDist(v5, v4);
        //          ▽C with respect to v5: (v5-v4) / dist(v5,v4)
        for(int currVertex = 0; currVertex < NUM_VERTEX; currVertex++){
            if(currVertex == 0 || currVertex == 6){
                continue;
            }
            for(int c = 0; c < C_PER_VERTEX; c++){
                // according to (18), each time, one element in the lambda vector is updated
                // For the toy constraint, ▽C is a sparse vector; for plant, it's the J
                int j = currVertex * C_PER_VERTEX + c;
                int relatedVertex = findRelateVertex(currVertex, c);
                if(relatedVertex == -1){
                    continue;
                }
                //compute delta_lambda using Eq(18)
                double eval_constraint = evalConstraint(currVertex, c);
                double scalar = evalScalar(currVertex, j, c);
                double alpha_tilde_j = 1 / (seconds * seconds) * alpha;
                double lambda_j = lambda(0,j);
                double delta_lambda = (-eval_constraint-alpha_tilde_j*lambda_j) / (scalar + alpha_tilde_j);

                // compute delta_x using Eq(17);
                //  since only one element of lambda is updated, so dx = M-1* ▽Cj(xi) * delta_lambda
                Eigen::VectorXd dxVector = evalGradient(currVertex, j, c) * delta_lambda;
                // convert VectorXd to MatrixXd as row-major
                Eigen::MatrixXd dx = Eigen::MatrixXd::Zero(3, NUM_VERTEX);
                int cnt = 0;
                for(int i = 0; i < NUM_VERTEX; i++){
                    for(int j = 0; j < 3; j++){
                        dx(j, i) = dxVector(cnt);
                        cnt++;
                    }
                }

                //update lambda
                lambda(j) += delta_lambda;
                //update x_tilde
                x_tilde += dx;
            }
        }

    }
    x_tilde.col(0) = x0.col(0);
    x_tilde.col(6) = x0.col(6);
    Eigen::MatrixXd v_next = 1 / seconds * (x_tilde - x);
    x = x_tilde;
    v = v_next;
    for(int i = 0; i < NUM_VERTEX; i++){
        vertices[i] = x.col(i);
    }
    m_shape.setVertices(vertices);
}

void Simulation::draw(Shader *shader)
{
    m_shape.draw(shader);
    m_ground.draw(shader);
}

void Simulation::toggleWire()
{
    m_shape.toggleWireframe();
}

void Simulation::initGround()
{
    std::vector<Vector3d> groundVerts;
    std::vector<Vector3i> groundFaces;
    groundVerts.emplace_back(-5, 0, -5);
    groundVerts.emplace_back(-5, 0, 5);
    groundVerts.emplace_back(5, 0, 5);
    groundVerts.emplace_back(5, 0, -5);
    groundFaces.emplace_back(0, 1, 2);
    groundFaces.emplace_back(0, 2, 3);
    m_ground.init(groundVerts, groundFaces);
}

int findRelateVertex(int currVertex, int c){
    // special logic for the toy cloth
    // 6  4  2  0
    // 7  5  3  1
    // 11 10 9  8
    // 15 14 13 12
    int up[16] = {-1, 0, -1, 2, -1, 4, -1, 6, 1, 3, 5, 7, 8, 9, 10, 11};
    int down[16] = {1, 8, 3, 9, 5, 10, 7, 11, 12 ,13, 14, 15, -1, -1, -1, -1};
    int left[16] = {2, 3, 4, 5, 6, 7, -1, -1, 9, 10, 11, -1, 13, 14, 15, -1};
    int right[16] = {-1, -1, 0, 1, 2, 3, 4, 5, -1, 8, 9, 10, -1, 12, 13, 14};
    int upleft[16] = {-1, 2, -1, 4, -1, 6, -1, -1, 3, 5, 7, -1, 9, 10, 11, -1};
    int upright[16] = {-1, -1, -1, 0, -1, 2, -1, 4, -1, 1, 3, 5, -1, 8, 9, 10};
    int downleft[16] = {3, 9, 5, 10, 7, 11, -1, -1, 13, 14, 15, -1, -1, -1, -1, -1};
    int downright[16] = {-1, -1, 1, 8, 3, 9, 5, 10, -1, 12, 13, 14, -1, -1, -1, -1};
    switch(c){
    case 0:
        return up[currVertex];
    case 1:
        return down[currVertex];
    case 2:
        return left[currVertex];
    case 3:
        return right[currVertex];
    case 4:
        return upleft[currVertex];
    case 5:
        return upright[currVertex];
    case 6:
        return downleft[currVertex];
    case 7:
        return downright[currVertex];
    default:
        return -1;
    }
}

double evalConstraint(int currVertex, int c){
    int relatedVertex = findRelateVertex(currVertex, c);
    // A basic distrance constraint
    Eigen::MatrixXd mm = x.col(currVertex);
    double curr_distance = (x.col(currVertex) - x.col(relatedVertex)).norm();
    double rest_distance = (x0.col(currVertex) - x0.col(relatedVertex)).norm();
    return curr_distance - rest_distance;
}

// scalar = ▽Cj M-1 ▽CjT
double evalScalar(int currVertex, int j, int c){
    Eigen::VectorXd gradient = evalGradient(currVertex, j, c);
    return gradient.transpose() * gradient;
}

Eigen::VectorXd evalGradient(int currVertex, int j, int c){
    Eigen::VectorXd gradient = Eigen::VectorXd::Zero(3 * NUM_VERTEX);
    // A basic gradient of distrance constraint for toy cloth
    int relatedVertex = findRelateVertex(currVertex, c);
    Eigen::VectorXd gradient_at_currVertex = x.col(currVertex) - x.col(relatedVertex);
    gradient_at_currVertex.normalize();
    gradient(3 * currVertex + 0) = gradient_at_currVertex(0);
    gradient(3 * currVertex + 1) = gradient_at_currVertex(1);
    gradient(3 * currVertex + 2) = gradient_at_currVertex(2);
    if(relatedVertex != 0 || relatedVertex != 6){
        gradient(3 * relatedVertex + 0) = -gradient_at_currVertex(0);
        gradient(3 * relatedVertex + 1) = -gradient_at_currVertex(1);
        gradient(3 * relatedVertex + 2) = -gradient_at_currVertex(2);
    }
    return gradient;
}
