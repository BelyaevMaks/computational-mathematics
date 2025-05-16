#include <GL/glut.h>
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <algorithm> 
#include <iomanip>   

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const double EPSILON = 1e-9;

struct State {
    glm::dvec3 position;
    glm::dvec3 momentum;
    glm::dquat quaternion;
    glm::dvec3 angular_momentum;
    double integrated_potential_energy_buoyancy; 

    State operator*(double scalar) const {
        State result;
        result.position = this->position * scalar;
        result.momentum = this->momentum * scalar;
        result.quaternion.x = this->quaternion.x * scalar;
        result.quaternion.y = this->quaternion.y * scalar;
        result.quaternion.z = this->quaternion.z * scalar;
        result.quaternion.w = this->quaternion.w * scalar;
        result.angular_momentum = this->angular_momentum * scalar;
        result.integrated_potential_energy_buoyancy = this->integrated_potential_energy_buoyancy * scalar;
        return result;
    }

    State operator+(const State& other) const {
        State result;
        result.position = this->position + other.position;
        result.momentum = this->momentum + other.momentum;
        result.quaternion = this->quaternion + other.quaternion;
        result.angular_momentum = this->angular_momentum + other.angular_momentum;
        result.integrated_potential_energy_buoyancy = this->integrated_potential_energy_buoyancy + other.integrated_potential_energy_buoyancy;
        return result;
    }
};

struct Constants {
    double mass;
    glm::dmat3 inertia_tensor; 
    double sphere_radius;      
};

struct ExternalForces {
    glm::dvec3 gravity;        
    glm::dvec3 drag_coefficient; 
    double fluid_density;      
};

State current_state;
Constants sim_constants;
ExternalForces ext_forces;
double simulation_time = 0.0;
double dt = 0.01; 

void initOpenGL();
void drawScene();
void reshape(int w, int h);
void update(int value);
State compute_derivatives(const State& state_arg, const Constants& consts_arg, const ExternalForces& forces_arg);
void runge_kutta_4(State& state_ref, const Constants& consts_ref, const ExternalForces& forces_ref, double time_step);
void print_energy(const State& state_arg, const Constants& consts_arg, const ExternalForces& forces_arg);
void midpoint_step(State& state, const Constants& constants, const ExternalForces& forces, double dt);

State compute_derivatives(const State& state_arg, const Constants& consts_arg, const ExternalForces& forces_arg) {
    State derivatives;

    glm::dvec3 velocity = state_arg.momentum / consts_arg.mass;
    derivatives.position = velocity;

    glm::dvec3 total_force = consts_arg.mass * forces_arg.gravity;

    double water_level = -4.0; 
    double sphere_bottom_y = state_arg.position.y - consts_arg.sphere_radius;
    double sphere_top_y = state_arg.position.y + consts_arg.sphere_radius;
    
    double submerged_volume = 0.0;
    glm::dvec3 buoyancy_force_vector = glm::dvec3(0.0);
    glm::dvec3 buoyancy_application_point_relative = glm::dvec3(0.0);

    if (sphere_top_y <= water_level) {
        submerged_volume = (4.0/3.0) * M_PI * consts_arg.sphere_radius * consts_arg.sphere_radius * consts_arg.sphere_radius;
        buoyancy_force_vector = glm::dvec3(0.0, forces_arg.fluid_density * 9.81 * submerged_volume, 0.0);
        buoyancy_application_point_relative = glm::dvec3(0.0); 
    } else if (sphere_bottom_y < water_level) {
        double h = water_level - sphere_bottom_y; 
        h = std::max(0.0, std::min(h, 2.0 * consts_arg.sphere_radius)); 

        if (h > EPSILON) {
            submerged_volume = (M_PI * h * h * (3.0 * consts_arg.sphere_radius - h)) / 3.0;
            buoyancy_force_vector = glm::dvec3(0.0, forces_arg.fluid_density * 9.81 * submerged_volume, 0.0);
            
            double y_centroid_from_base;
            if (std::abs(3.0 * consts_arg.sphere_radius - h) < EPSILON) { 
                 y_centroid_from_base = h / 2.0; 
            } else {
                 y_centroid_from_base = (h * (4.0 * consts_arg.sphere_radius - h)) / (4.0 * (3.0 * consts_arg.sphere_radius - h));
            }
            double absolute_buoyancy_center_y = sphere_bottom_y + y_centroid_from_base;
            buoyancy_application_point_relative = glm::dvec3(0.0, absolute_buoyancy_center_y - state_arg.position.y, 0.0);
        }
    }
    total_force += buoyancy_force_vector;
    glm::dvec3 drag_force = -forces_arg.drag_coefficient * velocity; 
    total_force += drag_force;
    derivatives.momentum = total_force;

    glm::dmat3 rotation_matrix = glm::toMat3(glm::normalize(state_arg.quaternion)); 
    glm::dmat3 inertia_tensor_world = rotation_matrix * consts_arg.inertia_tensor * glm::transpose(rotation_matrix);
    glm::dmat3 inv_inertia_tensor_world = glm::inverse(inertia_tensor_world);
    glm::dvec3 angular_velocity_world = inv_inertia_tensor_world * state_arg.angular_momentum;

    glm::dvec3 total_torque = glm::dvec3(0.0);
    glm::dvec3 r_buoyancy = buoyancy_application_point_relative; 
    total_torque += glm::cross(r_buoyancy, buoyancy_force_vector);
    derivatives.angular_momentum = total_torque;

    glm::dquat omega_quat(0.0, angular_velocity_world.x, angular_velocity_world.y, angular_velocity_world.z);
    derivatives.quaternion = 0.5 * state_arg.quaternion * omega_quat;
    /*
    Производная для интегрированной потенциальной энергии плавучести
    dU_A/dt = -F_A_y * v_y
    buoyancy_force_vector.y - это F_A_y (вверх)
    velocity.y - это v_y (derivatives.position.y)
    */
    derivatives.integrated_potential_energy_buoyancy = -buoyancy_force_vector.y * velocity.y;
    
    return derivatives;
}

void midpoint_step(State& state, const Constants& constants, const ExternalForces& forces, double dt) {
    // Вычисляем производные в начальной точке
    State k1 = compute_derivatives(state, constants, forces);
    
    // Находим состояние в средней точке (t + dt/2)
    State midpoint_state = state;
    midpoint_state.position += k1.position * (dt/2.0);
    midpoint_state.momentum += k1.momentum * (dt/2.0);
    midpoint_state.angular_momentum += k1.angular_momentum * (dt/2.0);
    
    // Обновление кватерниона в средней точке
    midpoint_state.quaternion.x += k1.quaternion.x * (dt/2.0);
    midpoint_state.quaternion.y += k1.quaternion.y * (dt/2.0);
    midpoint_state.quaternion.z += k1.quaternion.z * (dt/2.0);
    midpoint_state.quaternion.w += k1.quaternion.w * (dt/2.0);
    midpoint_state.quaternion = glm::normalize(midpoint_state.quaternion);
    
    midpoint_state.integrated_potential_energy_buoyancy += k1.integrated_potential_energy_buoyancy * (dt/2.0);
    
    // Вычисляем производные в средней точке
    State k2 = compute_derivatives(midpoint_state, constants, forces);
    
    // Обновляем состояние, используя производные в средней точке
    state.position += k2.position * dt;
    state.momentum += k2.momentum * dt;
    state.angular_momentum += k2.angular_momentum * dt;
    
    state.quaternion.x += k2.quaternion.x * dt;
    state.quaternion.y += k2.quaternion.y * dt;
    state.quaternion.z += k2.quaternion.z * dt;
    state.quaternion.w += k2.quaternion.w * dt;
    state.quaternion = glm::normalize(state.quaternion);
    
    state.integrated_potential_energy_buoyancy += k2.integrated_potential_energy_buoyancy * dt;
}
void print_energy(const State& state_arg, const Constants& consts_arg, const ExternalForces& forces_arg) {
    double translational_ke = 0.5 * consts_arg.mass * glm::dot(state_arg.momentum / consts_arg.mass, state_arg.momentum / consts_arg.mass);
    glm::dmat3 rotation_matrix = glm::toMat3(glm::normalize(state_arg.quaternion));
    glm::dmat3 inertia_tensor_world = rotation_matrix * consts_arg.inertia_tensor * glm::transpose(rotation_matrix);
    glm::dmat3 inv_inertia_tensor_world = glm::inverse(inertia_tensor_world);
    glm::dvec3 angular_velocity_world = inv_inertia_tensor_world * state_arg.angular_momentum;
    double rotational_ke = 0.5 * glm::dot(state_arg.angular_momentum, angular_velocity_world);
    double total_kinetic_energy = translational_ke + rotational_ke;

    double gravitational_potential_energy = consts_arg.mass * 9.81 * state_arg.position.y;
    
    double actual_buoyancy_potential_energy = state_arg.integrated_potential_energy_buoyancy;
    
    // Полная энергия: E = KE + PE_grav + U_A 
    double total_energy = total_kinetic_energy + gravitational_potential_energy + actual_buoyancy_potential_energy; 
    
    std::cout << "Time: " << simulation_time << std::fixed << std::setprecision(2)
              << ", KE: " << total_kinetic_energy 
              << ", PE_grav: " << gravitational_potential_energy 
              << ", PE_buoy_integ: " << actual_buoyancy_potential_energy // Измененный вывод
              << ", Total E: " << total_energy << std::endl;
}

void initOpenGL() {
    glClearColor(0.1f, 0.1f, 0.2f, 1.0f); 
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING); 
    glEnable(GL_LIGHT0);   

    GLfloat light_ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    GLfloat light_diffuse[] = { 0.8f, 0.8f, 0.8f, 1.0f };
    GLfloat light_specular[] = { 0.5f, 0.5f, 0.5f, 1.0f };
    GLfloat light_position[] = { 10.0f, 10.0f, 10.0f, 0.0f }; 

    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);

    glEnable(GL_COLOR_MATERIAL); 
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glShadeModel(GL_SMOOTH); 
}

void drawScene() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    gluLookAt(0.0, 2.0, 15.0,  
              0.0, 0.0, 0.0,   
              0.0, 1.0, 0.0);  

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_LIGHTING); 
    
    glColor4f(0.2f, 0.4f, 0.8f, 0.4f);  
    glPushMatrix();
    glTranslatef(0.0f, (float)(-4.0), 0.0f);  
    glBegin(GL_QUADS);
        glVertex3f(-15.0f, 0.0f, -15.0f);
        glVertex3f(-15.0f, 0.0f,  15.0f);
        glVertex3f( 15.0f, 0.0f,  15.0f);
        glVertex3f( 15.0f, 0.0f, -15.0f);
    glEnd();
    glPopMatrix();
    
    glEnable(GL_LIGHTING); 
    glDisable(GL_BLEND);

    glPushMatrix();
    glTranslatef(current_state.position.x, current_state.position.y, current_state.position.z);
    
    glm::dquat q_norm = glm::normalize(current_state.quaternion); 
    double angle_rad = 2.0 * acos(q_norm.w); 
    double angle_deg = glm::degrees(angle_rad);
    glm::dvec3 axis = glm::dvec3(q_norm.x, q_norm.y, q_norm.z);
    if (glm::length(axis) > EPSILON) { 
        axis = glm::normalize(axis);
        glRotatef(angle_deg, axis.x, axis.y, axis.z);
    }
    
    GLfloat sphere_ambient[] = {0.8f, 0.2f, 0.2f, 1.0f}; 
    GLfloat sphere_diffuse[] = {0.8f, 0.2f, 0.2f, 1.0f};
    GLfloat sphere_specular[] = {0.9f, 0.9f, 0.9f, 1.0f};
    GLfloat sphere_shininess[] = {50.0f};

    glMaterialfv(GL_FRONT, GL_AMBIENT, sphere_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, sphere_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, sphere_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, sphere_shininess);

    glutSolidSphere(sim_constants.sphere_radius, 30, 30); 
    glPopMatrix();
    
    glutSwapBuffers();
}
void update(int value) {
    midpoint_step(current_state, sim_constants, ext_forces, dt);
    simulation_time += dt;

    print_energy(current_state, sim_constants, ext_forces);
    glutPostRedisplay();
    glutTimerFunc(16, update, 0);
}

void reshape(int w, int h) {
    if (h == 0) h = 1; 
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, (double)w / (double)h, 0.1, 200.0); 
    glMatrixMode(GL_MODELVIEW);
}



int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1024, 768); 
    glutCreateWindow("Physics Simulation: Floating Sphere");

    // Инициализация состояния тела
    current_state.position = glm::dvec3(0.0, 0.0, 0.0); 
    current_state.momentum = glm::dvec3(250.0, 150.0, 0.0);  
    current_state.quaternion = glm::dquat(1.0, 0.0, 0.0, 0.0); 
    current_state.angular_momentum = glm::dvec3(0.0, 0.0, 15.0); 
    current_state.integrated_potential_energy_buoyancy = 0.0; // Инициализация новой переменной

    // Инициализация констант тела
    sim_constants.mass = 300.0; 
    sim_constants.sphere_radius = 0.5; 
    double I = (2.0 / 5.0) * sim_constants.mass * sim_constants.sphere_radius * sim_constants.sphere_radius;
    sim_constants.inertia_tensor = glm::dmat3(
        I, 0.0, 0.0, 
        0.0, I, 0.0, 
        0.0, 0.0, I  
    );

    // Инициализация внешних сил и параметров среды
    ext_forces.gravity = glm::dvec3(0.0, -9.81, 0.0); 
    ext_forces.fluid_density = 1000.0; 
    ext_forces.drag_coefficient = glm::dvec3(0.5, 0.5, 0.5); 


    std::cout << std::fixed << std::setprecision(2); 

    initOpenGL();
    glutDisplayFunc(drawScene);
    glutReshapeFunc(reshape);
    glutTimerFunc(0, update, 0); 
    glutMainLoop();
    return 0;
}

