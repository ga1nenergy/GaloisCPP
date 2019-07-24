//
// Created by ga1nenergy on 23.07.2019.
//

#include "remainder_code.h"
#include "gfalgorithms.h"

#include <algorithm>
#include <cmath>

using namespace galoiscpp;

namespace remainder_code {

    std::vector<std::vector<GFelement>> parity_check_matrix(const std::vector<GFelement> &locators, int t) {
        std::vector<std::vector<GFelement>> H;

        for (int i = 0; i < t; i++) {
            std::vector<GFelement> row = locators;
            for (auto & locator : row)
                locator = locator.inverse().power(i+1);
            H.push_back(row);
        }

        return H;
    }

    std::vector<std::vector<GFelement>> parity_check_matrix_punct(const std::vector<GFelement> &locators, scenario_t sc) {
        GaloisField *field = locators[0].getField();
        std::vector<std::vector<GFelement>> H;

        for (int i = 0; i < sc.t_strong; i++) {
            std::vector<GFelement> row(locators.size());
            if (i < sc.t_weak) {
                for (size_t j = 0; j < row.size(); j++) {
                    row[j] = locators[j].inverse().power(i + 1);
                }
            } else {
                for (size_t j = 0; j < row.size(); j++) {
                    if (j < sc.n_strong) {
                        row[j] = locators[j].inverse().power(i + 1);
                    } else {
                        row[j] = GFelement(field, 0);
                    }
                }
            }
            H.push_back(row);
        }

        return H;
    }

    GFpoly find_remainder(const std::vector<int> &msg, const std::vector<std::vector<GFelement>> &H) {
        if (msg.size() != H[0].size())
            throw std::logic_error("Incompatible sizes");

        std::vector<GFelement> syndrome(H.size());
        for (size_t i = 0; i < syndrome.size(); i++) {
            syndrome[i] = GFelement::dotint(H[i], msg);
        }

        GFpoly rem = newton_identitites(syndrome);
        return rem;
    }

    std::vector<GFpoly> find_strong_remainders(const std::vector<int> &msg, const std::vector<GFelement> &locators_strong,
                                               const std::vector<std::vector<GFelement>> &H,
                                               const std::vector<std::vector<GFelement>> &combs) {
        if (msg.size() != H[0].size())
            throw std::logic_error("Incompatible sizes");

        std::vector<GFelement> syndrome(H.size());
        for (size_t i = 0; i < syndrome.size(); i++) {
            syndrome[i] = GFelement::dotint(H[i], msg);
        }

        std::vector<GFpoly> strong_remainders;
        for (auto const & comb : combs) {
            std::copy(comb.begin(), comb.end(), syndrome.begin());

            for (auto & e : syndrome) {
                std::cout << e << " ";
            }
            std::cout << std::endl;

            GFpoly strong_remainder = newton_identitites(syndrome);
            strong_remainders.push_back(strong_remainder);
        }

        return strong_remainders;
    }

    std::vector<std::vector<GFelement>> find_locator_combinations(const std::vector<GFelement> &locators, int len) {
//        std::vector<std::vector<GFelement>> combs(pow(locators.size(), len));
//        for (auto & comb : combs) {
//            comb = std::vector<GFelement>(len);
//        }
        std::vector<std::vector<GFelement>> combs;

        std::vector<GFelement> counter(len);
        locator_combination_iter(locators, combs, counter, 0);

        return combs;
    }

    void locator_combination_iter(const std::vector<GFelement> &locators, std::vector<std::vector<GFelement>> &combs,
                                  std::vector<GFelement> &counter, int elem_number) {
        for (size_t i = 0; i < locators.size(); i++) {
            counter[elem_number] = locators[i];
            if (elem_number == counter.size() - 1) {
                combs.push_back(counter);
            } else {
                locator_combination_iter(locators, combs, counter, elem_number + 1);
            }
        }
    }

    /* TODO
     * 1) add rem and H field compatibility check
     */
    std::vector<int> decode(const std::vector<int> &received_msg, const std::vector<GFelement> &locators,
                            const std::vector<std::vector<GFelement>> &H, const GFpoly &rem) {
        GaloisField *field = rem.getField();
        int t = static_cast<int>(H.size());

        GFpoly received_rem = find_remainder(received_msg, H);
        GFpoly g_x(field, t + 1); g_x[t + 1] = GFelement(field, 1);

        std::cout << "Remainder: " << rem << std::endl;
        std::cout << "Inverse remainder: " << rem.inverse(g_x) << std::endl;
        std::cout << "Received remainder: " << received_rem << std::endl;
        std::cout << "received_rem / rem: " << GFpoly::deconvmod(received_rem, rem, g_x) << std::endl;
        GFpoly error_poly = GFpoly::euclid(GFpoly::deconvmod(received_rem, rem, g_x), g_x, t);

        std::cout << "Error poly: " << error_poly << std::endl;

        auto rts = error_poly.roots();
        auto decoded_msg = received_msg;
        for (auto const & rt : rts) {
            auto idx = std::find(locators.begin(), locators.end(), rt);
            if (idx != locators.end()) {
                decoded_msg[idx - locators.begin()]--;
            } else {
                std::cout << "Here throws Decoding error" << std::endl;
            }
        }

        return decoded_msg;
    }

    std::vector<std::vector<int>> decode_strong(const std::vector<int> &received_msg, const std::vector<GFelement> &locators,
                                  const std::vector<std::vector<GFelement>> &H, const GFpoly &rem,
                                  const std::vector<std::vector<GFelement>> &combs) {
        GaloisField *field = rem.getField();
        int t = static_cast<int>(H.size());

        auto strong_remainders = find_strong_remainders(received_msg, locators, H, combs);
        std::vector<std::vector<int>> decoded_messages;

        GFpoly g_x(field, t + 1); g_x[t + 1] = GFelement(field, 1);

        for (auto const & received_rem : strong_remainders) {
            std::cout << "Remainder: " << rem << std::endl;
            std::cout << "Inverse remainder: " << rem.inverse(g_x) << std::endl;
            std::cout << "Received remainder: " << received_rem << std::endl;
            std::cout << "received_rem / rem: " << GFpoly::deconvmod(received_rem, rem, g_x) << std::endl;
            GFpoly error_poly = GFpoly::euclid(GFpoly::deconvmod(received_rem, rem, g_x), g_x, t);

            std::cout << "Error poly: " << error_poly << std::endl;

            auto rts = error_poly.roots();
            if (rts.empty() && error_poly.getDegree() != 0) {
                std::cout << "Decoding error" << std::endl;
                continue;
            }
            auto decoded_msg = received_msg;
            int flag = 1;
            for (auto const & rt : rts) {
                auto idx = std::find(locators.begin(), locators.end(), rt);
                if (idx != locators.end()) {
                    decoded_msg[idx - locators.begin()]--;
                    if (decoded_msg[idx - locators.begin()] < 0) {
                        std::cout << "Decoding error" << std::endl;
                        flag = 0;
                        break;
                    }
                } else {
                    std::cout << "Decoding error" << std::endl;
                    flag = 0;
                    break;
                }
            }

            if (flag)
                decoded_messages.push_back(decoded_msg);
        }

        return decoded_messages;
    }

    std::vector<int> decode_multiuser(const std::vector<int> &received_msg, const std::vector<GFelement> &locators,
                                      const std::vector<std::vector<GFelement>> &H, const GFpoly &rem,
                                      const std::vector<std::vector<GFelement>> &combs, scenario_t sc) {
        std::vector<GFelement> locators_strong(locators.begin(), locators.begin() + sc.n_strong);

        std::cout << "Strong locators: ";
        for (auto & e : locators_strong) {
            std::cout << e << " ";
        }
        std::cout << std::endl;

        auto decoded_messages_strong = decode_strong(received_msg, locators_strong, H, rem, combs);

        for (auto & msg : decoded_messages_strong) {
            for (auto & e : msg) {
                std::cout << e << " ";
            }
            std::cout << std::endl;
        }

        std::vector<std::vector<GFelement>> H_weak(H.begin(), H.begin() + sc.t_weak);
        std::vector<std::vector<int>> all_decoded_messages;
        for (auto const & msg : decoded_messages_strong) {
            auto decoded_msg = decode(msg, locators, H_weak, rem);
            all_decoded_messages.push_back(decoded_msg);
        }

        int min_weight = 10000, min_idx = -1;
        for (size_t msg_idx = 0; msg_idx < all_decoded_messages.size(); msg_idx++) {
            int weight = 0;
            for (size_t i = 0; i < all_decoded_messages[msg_idx].size(); i++) {
                weight += abs(all_decoded_messages[msg_idx][i] - received_msg[i]);
            }

            if (weight < min_weight) {
                min_idx = msg_idx;
                min_weight = weight;
            }
        }

        return all_decoded_messages[min_idx];
    }

    std::vector<int> decode_multiuser_separable(const std::vector<int> &received_msg, const std::vector<GFelement> &locators,
                                                const std::vector<std::vector<GFelement>> &H,
                                                const std::vector<GFelement> &syndrome, scenario_t sc) {
        std::cout << "Syndrome: ";
        for (auto & e : syndrome) {
            std::cout << e << " ";
        }
        std::cout << std::endl;

        std::vector<GFelement> locators_strong(locators.begin(), locators.begin() + sc.n_strong);

        std::cout << "Strong locators: ";
        for (auto & e : locators_strong) {
            std::cout << e << " ";
        }
        std::cout << std::endl;

        auto H_strong = H;
        for (auto & row : H_strong) {
            row.resize(sc.n_strong);
        }

        std::cout << "Strong syndrome: ";
        auto strong_syndrome = syndrome;
        for (int i = 0; i < sc.t_weak; i++) {
            for (int j = sc.n_strong - 1; j < sc.n_all; j++)
                strong_syndrome[i] = strong_syndrome[i] - H[i][i];
            std::cout << strong_syndrome[i] << " ";
        }
        std::cout << std::endl;

        auto rem_strong = newton_identitites(strong_syndrome);
        std::cout << "Strong rem: " << rem_strong << std::endl;
        auto decoded_msg_strong = decode(received_msg, locators_strong, H_strong, rem_strong);

        std::vector<std::vector<GFelement>> H_weak(H.begin(), H.begin() + sc.t_weak - 1);
        std::vector<GFelement> syndrome_weak(syndrome.begin(), syndrome.begin() + sc.t_weak - 1);
        auto rem_weak = newton_identitites(syndrome_weak);

        std::cout << "Weak syndrome: ";
        for (auto & e : syndrome_weak) {
            std::cout << e << " ";
        }
        std::cout << std::endl;

        std::cout << "Weak H: ";
        for (auto & row : H_weak) {
            for (auto & e : row)
                std::cout << e << " ";
            std::cout << std::endl;
        }

        std::cout << "Weak rem: " << rem_weak << std::endl;

        auto decoded_msg = decode(decoded_msg_strong, locators, H_weak, rem_weak);

        return decoded_msg;
    }
};