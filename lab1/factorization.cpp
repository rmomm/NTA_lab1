#include <iostream>
#include "factorization.h"
#include "primality.h"
using namespace std;

uint64_t trialDivision(uint64_t n) {

    vector<uint64_t> primes = { 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47 };

    for (uint64_t p : primes) {
        if (n % p == 0)
        {
            cout << "Found divisor: " << p;
            return p;
        }
    }

    return 0;
}

uint64_t f(uint64_t x, uint64_t n) {
    return (mul_mod(x, x, n) + 1) % n;
}

uint64_t pollardRho(uint64_t n) {
    if (n % 2 == 0) {
        return 2;
    }

    uint64_t x = 2;
    uint64_t y = 2;
    uint64_t d = 1;

    while (d == 1) {
        x = f(x, n);
        y = f(f(y, n), n);

        uint64_t diff = (x > y) ? (x - y) : (y - x);
        d = gcd(diff, n);

        if (d == n) {
            return 0;
        }
    }

    return d;
}

int legendre(uint64_t a, uint64_t p) {
    uint64_t ls = pow_mod(a, (p - 1) / 2, p);
    if (ls == p - 1) {
        return -1;
    }
    return (int)ls;
}

vector<uint64_t> buildFactorBase(uint64_t n, uint64_t bound) {
    vector<uint64_t> base;
    base.push_back(2);
    for (uint64_t p = 3; p <= bound; p += 2) {
        bool prime = true;
        for (uint64_t d = 3; d * d <= p; d += 2) {
            if (p % d == 0) {
                prime = false;
                break;
            }
        }
        if (prime && legendre(n, p) == 1) {
            base.push_back(p);
        }
    }
    return base;
}

bool factorOverBase(uint64_t x, const vector<uint64_t>& base, vector<int>& exps) {
    exps.assign(base.size(), 0);
    for (size_t i = 0; i < base.size(); i++) {
        while (x % base[i] == 0) {
            x /= base[i];
            exps[i]++;
        }
    }
    return (x == 1);
}


vector<int> gaussianEliminationF2(vector<vector<int>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();


    vector<vector<int>> comb(rows, vector<int>(rows, 0));
    for (int i = 0; i < rows; i++) {
        comb[i][i] = 1;
    }

    int r = 0;
    for (int c = 0; c < cols && r < rows; c++) {
        int pivot = -1;
        for (int i = r; i < rows; i++) {
            if (matrix[i][c]) {
                pivot = i;
                break;
            }
        }
        if (pivot == -1) {
            continue;
        }


        swap(matrix[r], matrix[pivot]);
        swap(comb[r], comb[pivot]);


        for (int i = 0; i < rows; i++) {
            if (i != r && matrix[i][c]) {
                for (int j = c; j < cols; j++) {
                    matrix[i][j] ^= matrix[r][j];
                }
                for (int j = 0; j < rows; j++) {
                    comb[i][j] ^= comb[r][j];
                }
            }
        }
        r++;
    }


    for (int i = 0; i < rows; i++) {
        bool zero = true;
        for (int j = 0; j < cols; j++) {
            if (matrix[i][j]) {
                zero = false;
                break;
            }
        }
        if (zero) {

            vector<int> dep;
            for (int j = 0; j < rows; j++) {
                if (comb[i][j]) {
                    dep.push_back(j);
                }
            }
            return dep;
        }
    }
    return {};
}

uint64_t brillhartMorrison(uint64_t n) {
    uint64_t bound = 50000;
    auto base = buildFactorBase(n, bound);

    vector<vector<int>> matrix;
    vector<vector<int>> exponents;
    vector<uint64_t> bValues;

    uint64_t sqrt_n = (uint64_t)sqrt((double)n);
    if (sqrt_n * sqrt_n == n) { 
        return sqrt_n;
    }

    uint64_t v = 1;
    uint64_t u = sqrt_n;
    uint64_t a = sqrt_n;

    uint64_t b_pr = 1;
    uint64_t b = a % n;

    for (size_t step = 1; step < 1000000; step++) {
        uint64_t v_next = (n - u * u) / v;
        if (v_next == 0) { 
            break; 
        }

        int64_t Q_i = (step % 2 == 1) ? -(int64_t)v_next : (int64_t)v_next;

        vector<int> exps;
        if (factorOverBase((uint64_t)llabs(Q_i), base, exps)) {
            matrix.push_back(vector<int>(exps.size()));
            for (size_t j = 0; j < exps.size(); j++) {
                matrix.back()[j] = exps[j] % 2;
            }
            exponents.push_back(exps);
            bValues.push_back(b);
            if (matrix.size() >= base.size() + 50) { 
                break; 
            }
        }

        uint64_t a_next = (sqrt_n + u) / v_next;
        uint64_t u_next = a_next * v_next - u;
        uint64_t b_next = (mul_mod(a_next, b, n) + b_pr) % n;

        v = v_next;
        u = u_next;
        a = a_next;
        b_pr = b;
        b = b_next;
    }

    if (matrix.empty()) { 
        return 0; 
    }

    auto dep = gaussianEliminationF2(matrix);
    if (!dep.empty()) {
        uint64_t x = 1;
        vector<int> totalExp(base.size(), 0);

        for (int idx : dep) {
            x = mul_mod(x, bValues[idx], n);
            for (size_t j = 0; j < base.size(); j++) {
                totalExp[j] += exponents[idx][j];
            }
        }

        uint64_t y = 1;
        for (size_t j = 0; j < base.size(); j++) {
            y = mul_mod(y, pow_mod(base[j], totalExp[j] / 2, n), n);
        }

        uint64_t g = gcd((x + y) % n, n);
        if (g != 1 && g != n) { 
            return g; 
        }
        g = gcd((x + n - y) % n, n);
        if (g != 1 && g != n) { 
            return g; 
        }
    }

    return 0;
}




void factor(uint64_t n, vector<uint64_t>& factors) {
    if (n == 1) return;

    if (isPrime(n, 5)) {
        factors.push_back(n);
        return;
    }

    uint64_t d = trialDivision(n);

    if (d == 0) {
        d = pollardRho(n);
    }
    if (d == 0) {
        d = brillhartMorrison(n);
    }

    if (d == 0 || d == n) {
        factors.push_back(n);
        return;
    }

    factor(d, factors);
    factor(n / d, factors);
}

void printCanonical(uint64_t n, vector<uint64_t>& factors) {
    sort(factors.begin(), factors.end());

    cout << n << " = ";
    for (size_t i = 0; i < factors.size();) {
        uint64_t p = factors[i];
        int cnt = 0;
        while (i < factors.size() && factors[i] == p) {
            cnt++;
            i++;
        }
        cout << p;
        if (cnt > 1) cout << "^" << cnt;
        if (i < factors.size()) cout << " * ";
    }
    cout << "\n";
}