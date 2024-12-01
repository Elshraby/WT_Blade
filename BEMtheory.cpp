#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <unordered_map>
#include <algorithm>
#include <iomanip>

const double PI = 3.14159265358979323846;

struct AirfoilData {
    std::vector<double> angle_of_attack;
    std::vector<double> lift_coefficient;
    std::vector<double> drag_coefficient;
};

struct BladeElement {
    double radius;
    double chord;
    double twist;
    std::string airfoil_type;
};

struct BEMResult {
    double chord;
    double twist;
    double axial_induction;
    double angular_induction;
    double local_power_coefficient;
};

// Function prototypes
std::unordered_map<std::string, AirfoilData> parse_airfoil_data(const std::string& filename);
double interpolate(const std::vector<double>& x, const std::vector<double>& y, double x_val);
BEMResult bem_iteration(const BladeElement& element, const AirfoilData& airfoil, double wind_speed, double rotor_speed, int num_blades, double rotor_radius);

// Main BEM optimization function
std::vector<BEMResult> optimize_blade(const std::vector<BladeElement>& blade, 
                                      const std::unordered_map<std::string, AirfoilData>& airfoil_data,
                                      double wind_speed, double rotor_speed, int num_blades, double rotor_radius) {
    std::vector<BEMResult> results;
    
    for (const auto& element : blade) {
        auto it = airfoil_data.find(element.airfoil_type);
        if (it == airfoil_data.end()) {
            throw std::runtime_error("Unknown airfoil type: " + element.airfoil_type);
        }
        
        BEMResult result = bem_iteration(element, it->second, wind_speed, rotor_speed, num_blades, rotor_radius);
        results.push_back(result);
    }
    
    return results;
}


std::unordered_map<std::string, AirfoilData> parse_airfoil_data(const std::string& filename) {
    std::unordered_map<std::string, AirfoilData> airfoil_data;
    std::ifstream file(filename);
    std::string line;
    std::string current_airfoil;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        if (line.find("AIRFOIL_TYPE:") == 0) {
            current_airfoil = line.substr(13);
            airfoil_data[current_airfoil] = AirfoilData();
        } else if (line != "ANGLE_OF_ATTACK, LIFT_COEFFICIENT, DRAG_COEFFICIENT") {
            std::istringstream iss(line);
            std::string angle_str, lift_str, drag_str;
            
            if (std::getline(iss, angle_str, ',') && 
                std::getline(iss, lift_str, ',') && 
                std::getline(iss, drag_str)) {
                double angle = std::stod(angle_str);
                double lift = std::stod(lift_str);
                double drag = std::stod(drag_str);
                
                airfoil_data[current_airfoil].angle_of_attack.push_back(angle);
                airfoil_data[current_airfoil].lift_coefficient.push_back(lift);
                airfoil_data[current_airfoil].drag_coefficient.push_back(drag);
            }
        }
    }

    return airfoil_data;
}

double get_max_cl_cd_angle(const AirfoilData& data) {
    double max_cl_cd = -std::numeric_limits<double>::infinity();
    double best_angle = 0.0;

    for (size_t i = 0; i < data.angle_of_attack.size(); ++i) {
        // Avoid division by zero
        if (data.drag_coefficient[i] != 0) {
            double cl_cd = data.lift_coefficient[i] / data.drag_coefficient[i];
            if (cl_cd > max_cl_cd) {
                max_cl_cd = cl_cd;
                best_angle = data.angle_of_attack[i];
            }
        }
    }

    std::cout << "Max Cl/Cd: " << max_cl_cd << " at angle: " << best_angle << " degrees" << std::endl;
    return best_angle;
}

struct BladeInitializer {
    std::vector<double> radius;
    std::vector<double> chord;
    std::vector<double> twist;
    std::vector<std::string> airfoil_type;
};

BladeInitializer initialize_blade_geometry(double rotor_radius, int num_elements, double tip_speed_ratio, 
                                           int num_blades, const std::vector<std::string>& airfoil_types) { // should add the & in the vector strings??>?
    BladeInitializer initializer;
    
    // Constants
    const double alpha_design = 5.0 * M_PI / 180.0;  // Design angle of attack (5 degrees)
    const double cl_design = 0.8;  // Design lift coefficient

    // Initialize vectors
    initializer.radius.resize(num_elements);
    initializer.chord.resize(num_elements);
    initializer.twist.resize(num_elements);
    initializer.airfoil_type.resize(num_elements);

    // Calculate initial geometry
    for (int i = 0; i < num_elements; i++) {
        double r = rotor_radius * (0.2 + 0.8 * i / (num_elements - 1));  // Start at 20% of radius
        double lambda_r = tip_speed_ratio * r / rotor_radius; // This's too small will make the cord so biggg. So, It worked lol
        
        initializer.radius[i] = r;
        
        // Optimal relative wind angle
        double phi = (2.0 / 3.0) * atan(1/lambda_r);
        double F = (2/M_PI) * acos(exp((-num_blades/2 * (1-(r/rotor_radius)))/(r/rotor_radius * sin(phi))));
        // Optimal chord
        // initializer.chord[i] = (8.0 * M_PI * r * sin(phi)) / (3.0 * num_blades * cl_design * lambda_r);
        initializer.chord[i] = ((8.0 * M_PI * F * r * sin(phi)) / (3.0 * num_blades * cl_design)) * ((cos(phi) - (lambda_r * sin(phi)))/(sin(phi) + (lambda_r * cos(phi))));
        
        // Optimal twist (in degrees)
        initializer.twist[i] = (phi - alpha_design) * 180.0 / M_PI;

        // Assign airfoil type
        if (i < airfoil_types.size()) {
            initializer.airfoil_type[i] = airfoil_types[i];
        } else {
            // Use the last specified airfoil type for any remaining elements
            initializer.airfoil_type[i] = airfoil_types.back();
        }
    }

    return initializer;
}

double interpolate(const std::vector<double>& x, const std::vector<double>& y, double x_val) {
    auto it = std::lower_bound(x.begin(), x.end(), x_val);
    if (it == x.end()) return y.back();
    if (it == x.begin()) return y.front();
    
    int index = std::distance(x.begin(), it) - 1;
    double x1 = x[index];
    double x2 = x[index + 1];
    double y1 = y[index];
    double y2 = y[index + 1];
    
    return y1 + (x_val - x1) * (y2 - y1) / (x2 - x1);
}

BEMResult bem_iteration(const BladeElement& element, const AirfoilData& airfoil, 
                        double wind_speed, double rotor_speed, int num_blades, double rotor_radius) {
    const int MAX_ITERATIONS = 1000;
    const double CONVERGENCE_TOLERANCE = 1e-4;
    const double AIR_DENSITY = 1.225;  // kg/m^3
    
    double omega = rotor_speed * 2 * PI / 60;  // Convert rpm to rad/s
    double r = element.radius;
    double lambda_r = omega * r / wind_speed;
    double lambda = omega * rotor_radius / wind_speed;  // Tip speed ratio
    
    // Initial guess for induction factors
    double a = 0.3;
    double a_prime = 0.0;
    double cl, cd;
    
    double optimal_chord = element.chord;
    double optimal_twist = element.twist;
    
    for (int i = 0; i < MAX_ITERATIONS; ++i) {
        double phi = std::atan2(wind_speed * (1 - a), omega * r * (1 + a_prime));
        double best_angle = get_max_cl_cd_angle(airfoil);
        double alpha = (optimal_twist - best_angle) * PI / 180 + phi;
        
        cl = interpolate(airfoil.angle_of_attack, airfoil.lift_coefficient, alpha * 180 / PI);
        cd = interpolate(airfoil.angle_of_attack, airfoil.drag_coefficient, alpha * 180 / PI);
        
        double cn = cl * std::cos(phi) + cd * std::sin(phi);
        double ct = cl * std::sin(phi) - cd * std::cos(phi);
        
        // Prandtl's tip loss factor
        double f = (num_blades / 2) * (rotor_radius - r) / (r * std::sin(phi));
        double F = (2 / PI) * std::acos(std::exp(-f));
        F = std::max(F, 0.01);  // Prevent division by zero
        
        double sigma = optimal_chord * num_blades / (2 * PI * r);
        
        double a_new = 1 / (4 * F * std::sin(phi) * std::sin(phi) / (sigma * cn) + 1);
        double a_prime_new = 1 / (4 * F * std::sin(phi) * std::cos(phi) / (sigma * ct) - 1);
        
        // Calculate new optimal chord and twist
        //optimal_chord = (8 * PI * r * F) / (num_blades * cl * (1 - std::cos(phi)));
        optimal_twist = phi * 180 / PI - best_angle;
        
        std::cout << "Iteration " << i << ": r = " << r << ", phi = " << phi * 180 / PI 
                  << ", alpha = " << alpha * 180 / PI << ", cl = " << cl << ", cd = " << cd 
                  << ", a = " << a << ", a' = " << a_prime << std::endl;
        
        if (std::abs(a - a_new) < CONVERGENCE_TOLERANCE && std::abs(a_prime - a_prime_new) < CONVERGENCE_TOLERANCE) {
            // Converged
            double cp_local = 4 * a_new * (1 - a_new) * (1 - a_new) * lambda_r * lambda_r * F;
            
            return {optimal_chord, optimal_twist, a_new, a_prime_new, cp_local};
        }
        
        a = 0.75 * a + 0.25 * a_new;  // Relaxation factor to improve convergence
        a_prime = 0.75 * a_prime + 0.25 * a_prime_new;
    }
    
    // If we reach here, we didn't converge
    std::cerr << "Warning: BEM iteration did not converge for r = " << r << std::endl;
    
    // Return the last calculated values
    double cp_local = 4 * a * (1 - a) * (1 - a) * lambda_r * lambda_r * (r / rotor_radius) * (r / rotor_radius);
    return {optimal_chord, optimal_twist, a, a_prime, cp_local};
}

// Helper function to transform the blade
template<typename F>
std::vector<BladeElement> transform_blade(const std::vector<BladeElement>& blade, F transform) {
    std::vector<BladeElement> new_blade;
    std::transform(blade.begin(), blade.end(), std::back_inserter(new_blade), transform);
    return new_blade;
}

// Helper function to calculate total Cp
double calculate_total_cp(const std::vector<BladeElement>& blade, const std::vector<BEMResult>& results, double rotor_radius) {
    double total_cp = 0.0;
    for (size_t i = 0; i < blade.size() - 1; ++i) {
        double dr = blade[i+1].radius - blade[i].radius;
        double r = (blade[i+1].radius + blade[i].radius) / 2;
        total_cp += results[i].local_power_coefficient * 2 * r * dr / (rotor_radius * rotor_radius);
    }
    return total_cp;
}

int main() {
    try {
        // Parse airfoil data
        auto airfoil_data = parse_airfoil_data("airfoil_data.txt");
        
        // Set turbine parameters
        double wind_speed = 12.0;  // m/s
        double rotor_speed = 625;  // rpm
        int num_blades = 3;
        double rotor_radius = 0.85;  // m
        int num_elements = 12;
        double tip_speed_ratio = (rotor_speed * M_PI / 30) * rotor_radius / wind_speed;

        // Define airfoil types for different sections
        std::vector<std::string> airfoil_types = {
            "NACA4412", "NACA4412", "NACA4412", "NACA4412", // Root section
            "NACA2412", "NACA2412", "NACA2412", "NACA2412", // Mid section
            "NACA63415", "NACA63415", "NACA63415", "NACA63415", // Tip section
        };

        // Initialize blade geometry
        auto initializer = initialize_blade_geometry(rotor_radius, num_elements, tip_speed_ratio, num_blades, airfoil_types);

        // Create blade elements from the initializer
        std::vector<BladeElement> blade; 
        for (size_t i = 0; i < initializer.radius.size(); ++i) {
            blade.push_back({
                initializer.radius[i],
                initializer.chord[i],
                initializer.twist[i],
                initializer.airfoil_type[i]
            });
            // Debug print
            std::cout << "Blade Element " << i << ": Airfoil = " << initializer.airfoil_type[i] << std::endl;
        }
        
        // Optimize blade
        auto results = optimize_blade(blade, airfoil_data, wind_speed, rotor_speed, num_blades, rotor_radius);
        
        // Print results
        std::cout << "Radius(m)\tOptimal Chord(m)\tInitial Twist(deg)\tOptimal Twist(deg)\tAxial Ind.\tAngular Ind.\tLocal Cp\n";
        for (size_t i = 0; i < blade.size(); ++i) {
            std::cout << std::fixed << std::setprecision(3)
                      << blade[i].radius << "\t\t"
                      << results[i].chord << "\t\t\t"
                      << blade[i].twist << "\t\t\t"
                      << results[i].twist << "\t\t\t"
                      << results[i].axial_induction << "\t\t"
                      << results[i].angular_induction << "\t\t"
                      << results[i].local_power_coefficient << "\n";
        }
        
        // Calculate and print total power coefficient
        double total_cp = calculate_total_cp(blade, results, rotor_radius);
        std::cout << "\nTotal Power Coefficient: " << total_cp << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}