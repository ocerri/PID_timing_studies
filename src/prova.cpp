#include <iostream>
#include <vector>

using namespace std;

int main()
{
    // vector<int> v1 = {5, 6};
    // vector<int> v2;
    // cout << v1.empty() << endl;
    // cout << v2.empty() << endl;

    int a, b;
    cout << "Number?" << endl;
    cin >> a;
    cin >> b;
    cout << endl;

    double* arr = new double[a][b];

    for(unsigned int i = 0; i < a; i++)
    {
      for(unsigned int j = 0; i < b; i++)
      {
        cout << arr[i][j] << endl;
      }
    }

    return 0;
}
