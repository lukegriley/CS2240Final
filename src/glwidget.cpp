#include "glwidget.h"

#include <QApplication>
#include <QKeyEvent>
#include <iostream>

#include "plant/loader.h"

#define SPEED 1.5
#define ROTATE_SPEED 0.0025

using namespace std;

GLWidget::GLWidget(QWidget *parent) :
    QOpenGLWidget(parent),
    m_deltaTimeProvider(),
    m_intervalTimer(),
    m_sim(),
    m_camera(),
    m_shader(),
    m_forward(),
    m_sideways(),
    m_vertical(),
    m_lastX(),
    m_lastY(),
    m_capture(false)
{
    // GLWidget needs all mouse move events, not just mouse drag events
    setMouseTracking(true);

    // Hide the cursor since this is a fullscreen app
    QApplication::setOverrideCursor(Qt::ArrowCursor);

    // GLWidget needs keyboard focus
    setFocusPolicy(Qt::StrongFocus);

    // Function tick() will be called once per interva
    connect(&m_intervalTimer, SIGNAL(timeout()), this, SLOT(tick()));
}

GLWidget::~GLWidget()
{
    if (m_shader != nullptr) delete m_shader;
}

// ================== Basic OpenGL Overrides

void GLWidget::initializeGL()
{
    // Initialize GL extension wrangler
    glewExperimental = GL_TRUE;
    GLenum err = glewInit();
    if (err != GLEW_OK) fprintf(stderr, "Error while initializing GLEW: %s\n", glewGetErrorString(err));
    fprintf(stdout, "Successfully initialized GLEW %s\n", glewGetString(GLEW_VERSION));

    // Set clear color to white
    glClearColor(1, 1, 1, 1);

    // Enable depth-testing and backface culling
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    // Initialize the shader and simulation
    m_shader = new Shader(":/resources/shaders/shader.vert", ":/resources/shaders/shader.frag");
    m_plantShader = new Shader(":/resources/shaders/shader.vert", ":/resources/shaders/plant.frag");
    m_sim.init();

    // Initialize plant and renderer
    m_plant = load_plant("./data/plants/plant000.txt");
    m_plantRenderer.init(m_plant);

    // Initialize camera with a reasonable transform
    Eigen::Vector3f eye    = {0, 2, -5};
    Eigen::Vector3f target = {0, 1,  0};
    m_camera.lookAt(eye, target);
    m_camera.setOrbitPoint(target);
    m_camera.setPerspective(120, width() / static_cast<float>(height()), 0.1, 50);
    m_sim.loadCamera(&m_camera);

    m_deltaTimeProvider.start();
    m_intervalTimer.start(1000 / 120);
}

void GLWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    m_shader->bind();
    m_shader->setUniform("proj", m_camera.getProjection());
    m_shader->setUniform("view", m_camera.getView());
    m_sim.draw(m_shader);
    m_shader->unbind();

    m_plantShader->bind();
    m_plantShader->setUniform("proj", m_camera.getProjection());
    m_plantShader->setUniform("view", m_camera.getView());
    m_plantRenderer.render(m_plantShader);
    m_plantShader->unbind();
}

void GLWidget::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);
    m_camera.setAspect(static_cast<float>(w) / h);
}

// ================== Event Listeners

void GLWidget::mousePressEvent(QMouseEvent *event)
{
    m_capture = true;
    m_lastX = event->position().x();
    m_lastY = event->position().y();
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    if (!m_capture) return;

    int currX  = event->position().x();
    int currY  = event->position().y();

    int deltaX = currX - m_lastX;
    int deltaY = currY - m_lastY;

    if (deltaX == 0 && deltaY == 0) return;

    m_camera.rotate(deltaY * ROTATE_SPEED,
                    -deltaX * ROTATE_SPEED);

    m_lastX = currX;
    m_lastY = currY;
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
    m_capture = false;
}

void GLWidget::wheelEvent(QWheelEvent *event)
{
    float zoom = 1 - event->pixelDelta().y() * 0.1f / 120.f;
    m_camera.zoom(zoom);
}

void GLWidget::keyPressEvent(QKeyEvent *event)
{
    if (event->isAutoRepeat()) return;

    switch (event->key())
    {
    case Qt::Key_W: m_forward  += SPEED; break;
    case Qt::Key_S: m_forward  -= SPEED; break;
    case Qt::Key_A: m_sideways -= SPEED; break;
    case Qt::Key_D: m_sideways += SPEED; break;
    case Qt::Key_F: m_vertical -= SPEED; break;
    case Qt::Key_R: m_vertical += SPEED; break;
    case Qt::Key_C: m_camera.toggleIsOrbiting(); break;
    case Qt::Key_T: m_sim.toggleWire(); break;
    case Qt::Key_Z: m_sim.push(); break;
    case Qt::Key_X: m_sim.pull(); break;
    case Qt::Key_O: m_sim.toggleRunning(); break;
    case Qt::Key_P: m_sim.forwardOnce(); break;
    case Qt::Key_B: m_sim.goToBreakpoint(); break;
    case Qt::Key_Escape: QApplication::quit();
    }
}

void GLWidget::keyReleaseEvent(QKeyEvent *event)
{
    if (event->isAutoRepeat()) return;

    switch (event->key())
    {
    case Qt::Key_W: m_forward  -= SPEED; break;
    case Qt::Key_S: m_forward  += SPEED; break;
    case Qt::Key_A: m_sideways += SPEED; break;
    case Qt::Key_D: m_sideways -= SPEED; break;
    case Qt::Key_F: m_vertical += SPEED; break;
    case Qt::Key_R: m_vertical -= SPEED; break;
    }
}

// ================== Physics Tick

void GLWidget::tick()
{
    float deltaSeconds = m_deltaTimeProvider.restart() / 1000.f;
    double curr =0;
    while (curr < deltaSeconds) {
        double step = min(.003,deltaSeconds - curr);
        m_sim.update(step);
        curr += step;
    }

    m_plant.updateDiffusionDelta(deltaSeconds);
    this->m_plantRenderer.update_colors(m_plant);

    // Move camera
    auto look = m_camera.getLook();
    look.y() = 0;
    look.normalize();
    Eigen::Vector3f perp(-look.z(), 0, look.x());
    Eigen::Vector3f moveVec = m_forward * look.normalized() + m_sideways * perp.normalized() + m_vertical * Eigen::Vector3f::UnitY();
    moveVec *= deltaSeconds;
    m_camera.move(moveVec);
    Vector4f o_h = m_camera.getView() * Vector4f(0,0,0,1);
    Vector3f o(o_h[0],o_h[1],o_h[2]);
    Vector4f d_h = m_camera.getView() * Vector4f(m_camera.getLook()[0],m_camera.getLook()[1],m_camera.getLook()[2],0.f);
    Vector3f d(d_h[0],d_h[1],d_h[2]);
    m_sim.m_system.updateCameraPos(o,d);
    // Flag this view for repainting (Qt will call paintGL() soon after)
    update();
}
