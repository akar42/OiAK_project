#include <iostream>
#include <vector>
#include <tuple>
#include <numeric>
#include <random>
#include <cstdint>

using namespace std;

uint64_t convert_to_Montgomery_form(uint64_t number, uint64_t Gamma, uint64_t N)
{
	return (number * Gamma) % N;
}

uint64_t convert_from_Montgomery_form(uint64_t X_G, uint64_t Gamma_inv, uint64_t N)
{
	return (X_G * Gamma_inv) % N;
}

vector<int> RNS_representation(uint64_t X, const vector<int>& base)
{
	vector<int> remainders;

	for (int m: base)
	{
		remainders.push_back(X % m);
	}

	return remainders;
}

vector<int> MMM_in_RNS(const vector<int> X, const vector<int> Y, int N_hat, const int Gamma_inv, const vector<int> base)
{
	// Step 1
	vector<int> w;
	for (int i = 0; i < X.size(); i++)
	{
		w.push_back((X[i] * Y[i]) % base[i]);
	}

	// Step 2
	vector<int> N_hat_RNS = RNS_representation(N_hat, base);
	vector<int> q;
	for (int i = 0; i < X.size(); i++)
	{
		q.push_back((w[i] * N_hat_RNS[i]) % base[i]);
	}

	// Step 3
	int N = 1;
	for (int m: base)
	{
		N *= m;
	}

	vector<int> u;
	for (int i = 0; i < X.size(); i++)
	{
		u.push_back((q[i] * (N % base[i])) % base[i]);
	}

	// Step 4
	vector<int> v;
	for (int i = 0; i < X.size(); i++)
	{
		v.push_back((w[i] + u[i]) % base[i]);
	}

	// Step 5
	vector<int> Gamma_inv_RNS = RNS_representation(Gamma_inv, base);
	vector<int> Z;
	for (int i = 0; i < X.size(); i++)
	{
		Z.push_back((v[i] * Gamma_inv_RNS[i]) % base[i]);
	}

	return Z;
}


uint64_t multiplicative_inverse(uint64_t a, uint64_t m)
{
	int64_t m0 = m, t, q;
	int64_t x0 = 0, x1 = 1;

	if (m == 1)
		return 0;

	while (a > 1)
	{
		q = a / m;
		t = m;

		m = a % m;
		a = t;
		t = x0;

		x0 = x1 - q * x0;
		x1 = t;
	}

	if (x1 < 0)
		x1 += m0;

	return x1;
}

uint64_t chinese_remainder_theorem(const vector<int> &base, const vector<int> &remainders)
{
	uint64_t M = 1;
	for (int m : base)
	{
		M *= m;
	}

	uint64_t result = 0;

	for (int i = 0; i < base.size(); ++i)
	{
		uint64_t M_i = M / base[i];
		result += remainders[i] * multiplicative_inverse(M_i, base[i]) * M_i;
	}

	return result % M;
}

uint32_t generate_random_u_int_28_bit() {
    random_device rd;  
    mt19937 gen(rd()); 
    uniform_int_distribution<uint32_t> dis(0x08000000, 0x0FFFFFFF);

    return dis(gen);
}


int main()
{
	vector<int> base = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};

	uint64_t N = 1;

	for (int m: base)
	{
		N *= m;
	}

	cout << N << '\n';

	uint64_t Gamma = 7291011716708830237;

	uint64_t N_hat = multiplicative_inverse((-N + Gamma), Gamma);
	uint64_t Gamma_inv = multiplicative_inverse(Gamma, N);

	uint32_t X = generate_random_u_int_28_bit();
	uint32_t Y = generate_random_u_int_28_bit();

	cout << "X = " << X << '\t' << "Y =" << Y << '\n';

	uint64_t X_G = convert_to_Montgomery_form(X, Gamma, N);
	uint64_t Y_G = convert_to_Montgomery_form(Y, Gamma, N);

	vector<int> X_G_RNS = RNS_representation(X_G, base);
	vector<int> Y_G_RNS = RNS_representation(Y_G, base);

	vector<int> Z_G_RNS = MMM_in_RNS(X_G_RNS, Y_G_RNS, N_hat, Gamma_inv, base);

	uint64_t Z_G = chinese_remainder_theorem(base, Z_G_RNS);

	uint64_t Z = convert_from_Montgomery_form(Z_G, Gamma_inv, N);

	cout << "Z = X * Y = " << Z << '\n';

	return 0;
}
