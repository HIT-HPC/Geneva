#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector> 
#include <stdio.h>
#include <string.h>
#include <fstream> 
#include <iostream>
#include <regex>

typedef int dtype;
using namespace std;

void csr_to_matrix(dtype *value,dtype *colindex,dtype *rowptr,int n,int a,dtype** & M){
   M=new int*[n];
   for(int i=0;i<n;i++)
      M[i]=new int[n];
   for(int i=0;i<n;i++)
       for(int j=0;j<n;j++)
           M[i][j]=0;
   for(int i=0;i<n;i++)
       for(int j=rowptr[i];j<rowptr[i+1];j++)
           M[i][colindex[j]]=value[j];
   return;
}


int matrix_to_csr(int n, dtype **M, dtype* &value, dtype* & rowptr, dtype* & colindex){
    int i, j;
    int a = 0;
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            if(M[i][j] != 0)
                a++;
    value = new dtype[a];
    colindex = new int[a];
    rowptr = new int[n+1];
    int k = 0;
    int l = 0;
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++){
            if(j == 0)
                rowptr[l++] = k;
            if(M[i][j] != 0){
                value[k] = M[i][j];
                colindex[k] = j;
                k++;
            }
        }
   rowptr[l] = a;
   return a;
}

vector<string> split(const string& input, const string& regex) {
  
  std::regex re(regex);
  std::sregex_token_iterator first {input.begin(), input.end(), re, -1}, last;
  return {first, last};
}

void generate_sparse_matrix(dtype** & m, int n){
    m = new int*[n]();
    for(int i = 0; i < n; i++)
        m[i] = new int[n]();
    ifstream myfile("./DD.edges.format"); 
    string temp; 
    if (!myfile.is_open()) { 
        cout << "open file false." << endl; 
    }
    string comma = ",";
    while(getline(myfile,temp)) {
        auto res_str = split(temp, comma);
        int a = atoi(res_str[0].c_str());
        int b = atoi(res_str[1].c_str());
        m[a][b] = 1;
    }
    myfile.close();
}

void print_matrix(dtype **m, int n){
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++) {
            cout << m[i][j] << ",";
            if(j == n - 1)
                cout << endl;
        }
    return;
}



int main(){
    srand(time(0));
    int n;
    double s;
    cout << "n:";
    cin >> n;
    dtype **mat = NULL;
    dtype *value = NULL;
    int *colindex = NULL;
    int *rowptr = NULL;
    generate_sparse_matrix(mat, n);

    int a = matrix_to_csr(n, mat, value, rowptr, colindex);
    // dtype **mat_recover = NULL;
    // csr_to_matrix(value, colindex, rowptr, n, a, mat_recover);
    // cout << "matrix and csr transformation test" << endl;
    // int error = 0;
    // for(int i = 0; i < n; i++)
    //     for(int j = 0; j < n; j++)
    //        if(mat[i][j] != mat_recover[i][j])
    //            error = 1;
    // if(error == 1)
    //     cout << "test error!" << endl;
    // else
    //     cout << "test right!" << endl;
    cout << "a: " << a << endl;
    ofstream outfile("./DD.csr", ios::app);
    outfile << "value:" << endl;
    int value_num = sizeof(value) / sizeof(value[0]);
    for (int i = 0; i < value_num; i++){
        outfile << value[i] << endl;
    }
    outfile << "colindex:" << endl;
    int colinedx_num = sizeof(colindex) / sizeof(colindex[0]);
    for (int i = 0; i < colinedx_num; i++){
        outfile << colindex[i] << endl;
    }
    outfile << "rowptr:" << endl;
    int rowptr_num = sizeof(rowptr) / sizeof(rowptr[0]);
    for (int i = 0; i < rowptr_num; i++){
        outfile << rowptr[i] << endl;
    }
    return 0;
}

