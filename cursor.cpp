// Import necessary libraries
#include <iostream>
#include <vector>
#include <cmath>

// Define function to calculate Gaussian probability density function
double gaussian_pdf(double x, double mean, double stdev) {
    double exponent = -0.5 * pow((x - mean) / stdev, 2);
    double coefficient = 1 / (stdev * sqrt(2 * M_PI));
    return coefficient * exp(exponent);
}

// Define function to initialize parameters for EM algorithm
void initialize_parameters(std::vector<double>& weights, std::vector<double>& means, std::vector<double>& stdevs, int num_clusters) {
    // Initialize weights to be uniform
    double uniform_weight = 1.0 / num_clusters;
    for (int i = 0; i < num_clusters; i++) {
        weights.push_back(uniform_weight);
    }
    // Initialize means to be equally spaced
    double min_val = -10.0;
    double max_val = 10.0;
    double interval = (max_val - min_val) / (num_clusters - 1);
    for (int i = 0; i < num_clusters; i++) {
        means.push_back(min_val + i * interval);
    }
    // Initialize stdevs to be 1.0
    for (int i = 0; i < num_clusters; i++) {
        stdevs.push_back(1.0);
    }
}

// Define function to perform EM algorithm
void perform_em_algorithm(std::vector<double>& data, std::vector<double>& weights, std::vector<double>& means, std::vector<double>& stdevs, int num_clusters, int num_iterations) {
    int num_data_points = data.size();
    int num_params = num_clusters * 3;
    std::vector<double> responsibilities(num_clusters);
    std::vector<double> params(num_params);
    for (int i = 0; i < num_iterations; i++) {
        // E-step: calculate responsibilities
        for (int j = 0; j < num_data_points; j++) {
            double sum = 0.0;
            for (int k = 0; k < num_clusters; k++) {
                double pdf = gaussian_pdf(data[j], means[k], stdevs[k]);
                sum += weights[k] * pdf;
            }
            for (int k = 0; k < num_clusters; k++) {
                double pdf = gaussian_pdf(data[j], means[k], stdevs[k]);
                responsibilities[k] = weights[k] * pdf / sum;
            }
            // M-step: update parameters
            for (int k = 0; k < num_clusters; k++) {
                double sum_responsibilities = 0.0;
                double sum_weighted_data = 0.0;
                double sum_weighted_squared_data = 0.0;
                for (int l = 0; l < num_data_points; l++) {
                    sum_responsibilities += responsibilities[l];
                    sum_weighted_data += responsibilities[l] * data[l];
                    sum_weighted_squared_data += responsibilities[l] * pow(data[l], 2);
                }
                weights[k] = sum_responsibilities / num_data_points;
                means[k] = sum_weighted_data / sum_responsibilities;
                stdevs[k] = sqrt(sum_weighted_squared_data / sum_responsibilities - pow(means[k], 2));
            }
        }
    }
}

// Example usage
int main() {
    // Generate example data
    std::vector<double> data;
    for (int i = 0; i < 1000; i++) {
        double val = (i % 2 == 0) ? 5.0 : -5.0;
        val += rand() / (double)RAND_MAX;
        data.push_back(val);
    }
    // Initialize parameters
    std::vector<double> weights;
    std::vector<double> means;
    std::vector<double> stdevs;
    int num_clusters = 2;
    initialize_parameters(weights, means, stdevs, num_clusters);
    // Perform EM algorithm
    int num_iterations = 10;
    perform_em_algorithm(data, weights, means, stdevs, num_clusters, num_iterations);
    // Print results
    for (int i = 0; i < num_clusters; i++) {
        std::cout << "Cluster " << i << ":\n";
        std::cout << "Weight: " << weights[i] << "\n";
        std::cout << "Mean: " << means[i] << "\n";
        std::cout << "Standard deviation: " << stdevs[i] << "\n";
    }
    return 0;
}