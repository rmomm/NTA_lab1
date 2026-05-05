#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <string>
using namespace std;

uint64_t trialDivision(uint64_t n);
uint64_t f(uint64_t x, uint64_t n);
uint64_t pollardRho(uint64_t n);

int legendre(uint64_t a, uint64_t p);
vector<uint64_t> buildFactorBase(uint64_t n, uint64_t bound);
vector<uint64_t> buildFactorBase(uint64_t n, uint64_t bound);
bool factorOverBase(uint64_t x, const vector<uint64_t>& base, vector<int>& exps);
vector<int> gaussianEliminationF2(vector<vector<int>>& M);
uint64_t brillhartMorrison(uint64_t n);

void factor(uint64_t n, vector<uint64_t>& factors);
void printCanonical(uint64_t n, vector<uint64_t>& factors);