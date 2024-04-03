#include <iostream>

using namespace std;
void RNS_transformation(int x, int *base, int base_length, int array[]);
void RNS_multiplication(int A[], int B[], int base_length, int output[], int *base);
int modInverse(int a, int m);
int chineseRemainderTheorem(int n[], int a[], int base_length);

int main()
{
  int modular_base[] = {3, 4, 5};
  int numberA = 8;
  int numberB = 6;
  int A[3];
  int B[3];
  int Z[3];
  int zxc;
  RNS_transformation(numberA, modular_base, 3, A);
  RNS_transformation(numberB, modular_base, 3, B);
  RNS_multiplication(A, B, 3, Z, modular_base);
  zxc = chineseRemainderTheorem(modular_base, Z, 3);
  cout << zxc << endl;
  return 0;
}

void RNS_transformation(int x, int *base, int base_length, int array[])
{
  for (int i = 0; i < base_length; i++)
  {
    int r = x % base[i];
    array[i] = r;
  }
}

void RNS_multiplication(int A[], int B[], int base_length, int output[], int *base)
{
  for (int i = 0; i < base_length; i++)
  {
    output[i] = (A[i] * B[i]) % base[i];
  }
}

int modInverse(int a, int m)
{
  int m0 = m, t, q;
  int x0 = 0, x1 = 1;

  if (m == 1)
    return 0;

  while (a > 1)
  {
    q = a / m;
    t = m;

    m = a % m, a = t;

    t = x0;
    x0 = x1 - q * x0;
    x1 = t;
  }

  if (x1 < 0)
    x1 += m0;

  return x1;
}

int chineseRemainderTheorem(int n[], int a[], int base_length)
{
  int prod = 1;
  int result = 0;

  for (int i = 0; i < base_length; i++)
  {
    prod *= n[i];
  }

  for (int i = 0; i < base_length; i++)
  {
    int pp = prod / n[i];
    int inv = modInverse(pp, n[i]);
    if (inv == -1)
    {
      cerr << "Modular inverse does not exist for " << pp << " modulo " << n[i] << endl;
      return -1;
    }
    result += a[i] * inv * pp;
  }

  result = (result % prod + prod) % prod;
  return result;
}
