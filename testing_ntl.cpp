#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

int main()
{
   ZZ a, b, c;

   cin >> a;
   cin >> b;
   c = (a+1)*(b+1);
   cout << c << "\n";
}
