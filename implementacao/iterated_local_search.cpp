// Autor : MATEUS NUNES CAMPOS
// Cartão : 00268613
// Disciplina : OTIMIZAÇÃO COMBINATÓRIA
// Professor : MARCUS RITT
// Data : 19/08/2024

#include <algorithm>
#include <fstream>
#include <future>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <chrono>
#include <unordered_set>
#include <random>

//Hashing function for quick incompatibilities lookup
struct pair_hash {
    template<class T1, class T2>
    std::size_t operator ()(const std::pair<T1, T2> &p) const {
        return std::hash<T1>()(p.first) ^ std::hash<T2>()(p.second);
    }
};

//Parameters
int NEIGHBORHOOD = 20;
int MAX_ITERATIONS = 5000;
int NO_IMPROVEMENT_THRESHOLD = 15;
int MAX_ITERATIONS_WITHOUT_IMPROVEMENT = -1; //Disabled by default - use when time becomes too big of a problem
unsigned int SEED;

// Random number generator
std::mt19937 rng;

//Global variables
int total_ingredients, total_incompatibilities, max_weight;
std::vector<int> tastes_vector;
std::vector<int> weights_vector;
std::unordered_set<std::pair<int, int>, pair_hash> incompatibilities_set;

//Given two ingredients, check if they are incompatible
bool are_compatible(int ingredient1, int ingredient2) {
    return incompatibilities_set.find({ingredient1, ingredient2}) == incompatibilities_set.end() &&
           incompatibilities_set.find({ingredient2, ingredient1}) == incompatibilities_set.end();
}

//Must be a file in the format defined by the professor
bool parse_file(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file " << filename << std::endl;
        return false;
    }

    std::cout << "[+] File successfully open. Parsing data." << std::endl;

    std::istringstream line_buffer;
    std::string line;
    std::getline(file, line);
    line_buffer.str(line);

    line_buffer >> total_ingredients >> total_incompatibilities >> max_weight;

    std::cout << "[#] Total ingredients: " << total_ingredients << std::endl;
    std::cout << "[#] Total incompatibilities: " << total_incompatibilities << std::endl;
    std::cout << "[#] Max weight: " << max_weight << std::endl;

    // Skip empty line
    std::getline(file, line);

    // Parse tastes
    while (std::getline(file, line)) {
        if (line.empty()) {
            break; // Stop reading if an empty line is found
        }

        line_buffer.clear();
        line_buffer.str(line);

        int number = 0;
        while (line_buffer >> number) {
            tastes_vector.push_back(number);
        }
    }

    if (total_ingredients != tastes_vector.size()) {
        std::cerr << "[-] Error parsing data. Expected " << total_ingredients << " tastes, got " << tastes_vector.size()
                << std::endl;
        return false;
    }

    // Parse weights
    while (std::getline(file, line)) {
        if (line.empty()) {
            break; // Stop reading if an empty line is found
        }

        line_buffer.clear();
        line_buffer.str(line);

        int number = 0;
        while (line_buffer >> number) {
            weights_vector.push_back(number);
        }
    }

    if (total_ingredients != weights_vector.size()) {
        std::cerr << "[-] Error parsing data. Expected " << total_ingredients << " weights, got " << weights_vector.
                size() << std::endl;
        return false;
    }

    // Parse incompatibilities
    while (std::getline(file, line)) {
        if (line.empty()) {
            break; // Stop reading if an empty line is found
        }

        line_buffer.clear();
        line_buffer.str(line);

        int ingredient1 = 0;
        int ingredient2 = 0;
        while (line_buffer >> ingredient1 >> ingredient2) {
            incompatibilities_set.emplace(ingredient1, ingredient2);
        }
    }

    if (total_incompatibilities != incompatibilities_set.size()) {
        std::cerr << "[-] Error parsing data. Expected " << total_incompatibilities << " incompatibilities, got " <<
                incompatibilities_set.size() << std::endl;
        return false;
    }

    std::cout << "[+] Data parsed successfully. Now proceeding to find the best value." << std::endl;
    return true;
}


//Generates a deterministic initial solution by calculating a score based on the ratio of taste/weight and sorting higher->lower.
//The initial solution will be the first n elements of this sorted vector that are compatible and doesn't exceed the weight limit.

std::vector<int> generate_initial_solution() {
    std::vector<int> initial_solution;
    std::vector<std::pair<double, int> > ingredients_score;
    ingredients_score.reserve(total_ingredients);

    for (int i = 0; i < total_ingredients; i++) {
        double score = static_cast<double>(tastes_vector[i]) / weights_vector[i];
        ingredients_score.emplace_back(score, i);
    }

    std::sort(ingredients_score.begin(), ingredients_score.end(),
              [](const std::pair<double, int> &a, const std::pair<double, int> &b) {
                  return a.first > b.first;
              });

    std::unordered_set<int> used_ingredients;
    int current_weight = 0;

    // Construct the initial solution
    for (const auto &score_index: ingredients_score) {
        int index = score_index.second;
        int weight = weights_vector[index];

        if (current_weight + weight > max_weight) {
            continue; // Go to next ingredient if adding this ingredient exceeds max weight
        }

        bool can_add = true;
        for (const auto &selected: used_ingredients) {
            if (incompatibilities_set.find({selected, index}) != incompatibilities_set.end() ||
                incompatibilities_set.find({index, selected}) != incompatibilities_set.end()) {
                can_add = false;
                break;
            }
        }

        if (can_add) {
            initial_solution.push_back(index);
            used_ingredients.insert(index);
            current_weight += weight;
        }
    }

    return initial_solution;
}

// Debugging function for showing a vector of ingredients info
void print_vector(const std::vector<int> &vec) {
    std::cout << "[";
    if (!vec.empty()) {
        std::cout << vec[0];
        for (size_t i = 1; i < vec.size(); ++i) {
            std::cout << ", " << vec[i];
        }
    }
    std::cout << "]" << std::endl;

    int total_weight = 0;
    int total_score = 0;
    for (int index: vec) {
        total_weight += weights_vector[index];
        total_score += tastes_vector[index];
    }

    std::cout << "Total weight: " << total_weight << std::endl;
    std::cout << "Total taste: " << total_score << std::endl;
}


//Given a solution, returns its total weight
int get_solution_weight(const std::vector<int> &solution) {
    int total_weight = 0;

    for (int k: solution) {
        total_weight += weights_vector[k];
    }

    return total_weight;
}

//Given a solution, returns its total taste
int get_solution_taste(const std::vector<int> &solution) {
    int total_taste = 0;

    for (int k: solution) {
        total_taste += tastes_vector[k];
    }

    return total_taste;
}

//Given an ingredient and a solution, returns TRUE if the ingredient is present in the solution, FALSE otherwise
bool is_ingredient_in_solution(std::vector<int> &solution, const int ingredient) {
    for (const int k: solution) {
        if (solution[k] == ingredient) {
            return true;
        }
    }
    return false;
}

//Given a solution weight and an ingredient weight, returns TRUE if adding this ingredient will exceed the MAX_WEIGHT. FALSE otherwise.
bool will_overflow(const int solution_weight, const int ingredient_weight) {
    return (solution_weight + ingredient_weight > max_weight);
}

bool will_be_compatible(std::vector<int> &solution, const int ingredient) {
    return std::all_of(solution.begin(), solution.end(),
                       [&](const int k) {
                           return are_compatible(ingredient, k);
                       });
}


//Perturbates a solution based on the iterations_without_improving parameter. Removes 1*perturbation_strength random element(s).
std::vector<int> perturbate(std::vector<int> &solution, int iterations_without_improving) {
    if (solution.empty()) {
        std::cout << "[-] Solution is empty." << std::endl;
        return {};
    }

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<int> dist(0, solution.size() - 1);
    int perturbation_strength = (iterations_without_improving >= NO_IMPROVEMENT_THRESHOLD ? 4 : 2);
    //Increases strength if no improvements are found

    //Some stronger perturbations based on certain "breakpoints"
    if (iterations_without_improving == 50) {
        perturbation_strength =  (solution.size() - 1)/4;
    }

    if (iterations_without_improving == 99) {
        perturbation_strength = (solution.size() - 1)/3;
    }

    if (iterations_without_improving == 150 || iterations_without_improving == 200 || iterations_without_improving == 250) {
        perturbation_strength = (solution.size() - 1)/2; // half
    }

    for (int i = 0; i < perturbation_strength; ++i) {
        if (solution.empty()) {
            break;
        }
        int to_be_removed = dist(rng);
        solution.erase(solution.begin() + to_be_removed);
    }

    return solution;
}

bool already_exists_in_solution(std::vector<int> &solution, int ingredient) {
    if (std::find(solution.begin(), solution.end(), ingredient) != solution.end()) {
        return true;
    } else return false;
}


//Given two solutions, return the solution with best total taste
std::vector<int> get_best_solution(const std::vector<int> &solution1, const std::vector<int> &solution2) {
    const int solution1_score = get_solution_taste(solution1);
    const int solution2_score = get_solution_taste(solution2);

    if (solution1_score > solution2_score) {
        return solution1;
    } else {
        return solution2;
    }
}

//Given a solution, ingredient that were already removed (we don't want to add the same ingredient we just removed) and an history of ingredients, returns the first viable ingredient found
int find_viable_ingredient(std::vector<int> &solution, std::vector<int> history, int removed_ingredient) {
    const int solution_weight = get_solution_weight(solution);

    for (int i = 0; i < total_ingredients; i++) {
        if (already_exists_in_solution(solution, i)) {
            continue;
        }

        if (std::find(history.begin(), history.end(), i) != history.end()) {
            continue;
        }

        if (i == removed_ingredient) {
            continue;
        }

        if (will_overflow(solution_weight, weights_vector[i])) {
            continue;
        }

        if (!will_be_compatible(solution, i)) {
            continue;
        }

        // If all checks passed, return the ingredient
        return i;
    }


    return -1;
}

std::vector<int> available_ingredients(std::vector<int> &solution, int removed_ingredient) {
    std::vector<int> available_ingredients = {};

    while (true) {
        int ingredient = find_viable_ingredient(solution, available_ingredients, removed_ingredient);

        if (ingredient == -1) {
            break;
        }
        available_ingredients.push_back(ingredient);
    }
    return available_ingredients;
}

//Updates the available ingredients list without having to check for each ingredient possible
std::vector<int> update_available_list(std::vector<int> &available_list, int last_added_ingredient,
                                       int solution_weight) {
    if (available_list.empty()) {
        return {};
    }

    std::vector<int> updated_available_list = available_list;
    updated_available_list.erase(std::find(updated_available_list.begin(), updated_available_list.end(),
                                           last_added_ingredient));
    std::vector<int> to_remove;

    for (int k: updated_available_list) {
        bool is_incompatible = !are_compatible(k, last_added_ingredient);
        bool exceeds_weight = solution_weight + weights_vector[k] > max_weight;

        if (is_incompatible || exceeds_weight) {
            to_remove.push_back(k);
        }
    }

    // Remove elements that need to be removed
    for (int k: to_remove) {
        updated_available_list.erase(std::remove(updated_available_list.begin(), updated_available_list.end(), k),
                                     updated_available_list.end());
    }
    return updated_available_list;
}

//Generates a neighborhood of the size NEIGHBORHOOD - Parallelized for multithreading, may have some overhead
std::vector<std::vector<int> > generate_neighborhood(std::vector<int> &solution) {
    std::vector<std::future<std::vector<std::vector<int> > > > futures;
    int num_threads = std::thread::hardware_concurrency();
    int chunk_size = NEIGHBORHOOD / num_threads;
    std::vector<std::vector<int> > neighborhood;

    for (int i = 0; i < num_threads; i++) {
        int start = i * chunk_size;
        int end = (i == num_threads - 1) ? NEIGHBORHOOD : start + chunk_size;

        futures.push_back(std::async(std::launch::async, [=, &solution]() -> std::vector<std::vector<int> > {
            std::vector<std::vector<int> > local_neighborhood;
            std::random_device dev;
            std::mt19937 rng(dev());

            for (int k = start; k < end; k++) {
                std::uniform_int_distribution<int> dist(0, solution.size() - 1);
                int to_be_removed = dist(rng);

                std::vector<int> new_solution = solution;

                new_solution.erase(new_solution.begin() + to_be_removed);


                std::vector<int> viable_options = available_ingredients(new_solution, solution[to_be_removed]);


                if (!viable_options.empty()) {
                    std::uniform_int_distribution<int> dis(0, viable_options.size() - 1);
                    int to_be_added = viable_options[dis(rng)];
                    new_solution.push_back(to_be_added);

                    viable_options = update_available_list(viable_options, to_be_added,
                                                           get_solution_weight(new_solution));


                    while (get_solution_weight(new_solution) <= max_weight && !viable_options.empty()) {
                        std::uniform_int_distribution<int> dis2(0, viable_options.size() - 1);
                        int addAt = dis2(rng);

                        new_solution.push_back(viable_options[addAt]);
                        viable_options = update_available_list(viable_options, viable_options[addAt],
                                                               get_solution_weight(new_solution));
                    }
                }

                local_neighborhood.push_back(new_solution);
            }

            return local_neighborhood;
        }));
    }

    for (auto &f: futures) {
        auto local_neighborhood = f.get();
        neighborhood.insert(neighborhood.end(), local_neighborhood.begin(), local_neighborhood.end());
    }

    return neighborhood;
}

//Simple local search function
std::vector<int> local_search(std::vector<int> &solution, int no_improvement_count) {
    int best_objective_value = get_solution_taste(solution);
    std::vector<int> best_solution = solution;

    std::vector<std::vector<int> > neighborhood = generate_neighborhood(solution);


    for (const auto &neighbor: neighborhood) {
        const int neighbor_value = get_solution_taste(neighbor);

        //May accept worse solutions on checkpoints (for escaping local maxima)
        if (no_improvement_count == 35 || no_improvement_count == 50 || no_improvement_count == 75 ||
            no_improvement_count == 100 || no_improvement_count == 125 || no_improvement_count == 150 ||
            no_improvement_count == 175 || no_improvement_count == 200) {
            best_solution = neighbor;
            break;
        }

        //Best neighbor accepted
        if (neighbor_value > best_objective_value) {
            best_solution = neighbor;
            best_objective_value = neighbor_value;
        }
    }

    return best_solution;
}


std::vector<int> ILS() {
    std::vector<int> initial_solution = generate_initial_solution();


    std::vector<int> current_solution = initial_solution;
    std::vector<int> best_solution = initial_solution;
    int objective_function = get_solution_taste(best_solution);
    int iterations = 0;
    int no_improvement_count = 0;

    print_vector(initial_solution);


    std::cout << "[+] Running. " << std::endl;
    while (iterations < MAX_ITERATIONS) {
        if (no_improvement_count == MAX_ITERATIONS_WITHOUT_IMPROVEMENT) {
            std::cerr << "\r[-] Max iterations without improvement reached. Stopping." << std::endl;
            break;
        }

        std::cout << "\r[#] Iteration number : " << iterations << "/" << MAX_ITERATIONS << std::flush;
        //perturb
        std::vector<int> new_solution = perturbate(current_solution, no_improvement_count);


        //local optimum
        std::vector<int> local_optimum = local_search(new_solution, no_improvement_count);

        //update
        if (get_solution_taste(local_optimum) > objective_function) {
            best_solution = local_optimum;
            objective_function = get_solution_taste(best_solution);
            no_improvement_count = 0;
        } else {
            no_improvement_count++;
        }

        if (no_improvement_count == 100 && MAX_ITERATIONS_WITHOUT_IMPROVEMENT == -1) {
            no_improvement_count = 0;
        }


        current_solution = local_optimum;

        iterations++;
    }
    return best_solution;
}


int main(const int argc, char *argv[]) {
    if (argc < 2 || argc > 6) {
        std::cerr << "Usage: ./iterated_local_search <filename> \n"
                << "Optional arguments: \n"
                << "<neighborhood size (default 20)> \n"
                << "<max iterations (default 5000)> \n"
                << "<iterations with no improvement threshold (default 15)> \n"
                << "<max iterations without improvement (disabled by default)>"
                << std::endl;
        return 1;
    }

    // Parse optional arguments if provided
    if (argc >= 3) {
        NEIGHBORHOOD = std::stoi(argv[2]);
    }

    if (argc >= 4) {
        MAX_ITERATIONS = std::stoi(argv[3]);
    }

    if (argc >= 5) {
        NO_IMPROVEMENT_THRESHOLD = std::stoi(argv[4]);
    }

    if (argc == 6) {
        MAX_ITERATIONS_WITHOUT_IMPROVEMENT = std::stoi(argv[5]);
    }

    // Parse the input file
    if (!parse_file(argv[1])) {
        std::cerr << "Failed to parse the file: " << argv[1] << std::endl;
        return 1;
    }

    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();

    // Run the Iterated Local Search algorithm
    std::vector<int> final_solution = ILS();

    // Stop timing
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;

    // Convert elapsed time to minutes and seconds
    int minutes = static_cast<int>(elapsed_seconds.count()) / 60;
    int seconds = static_cast<int>(elapsed_seconds.count()) % 60;

    // Output the final solution
    std::cout << "\nFinal solution: " << std::endl;
    print_vector(final_solution);

    // Output the elapsed time
    std::cout << "Time taken: " << minutes << " minutes and " << seconds << " seconds." << std::endl;

    return 0;
}
