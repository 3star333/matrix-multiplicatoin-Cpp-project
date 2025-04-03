#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <random>
#include <chrono>
#include <string>

using namespace std;
using namespace chrono;
//functon prototypes:
//no classes because i hate myself and other people sorry (program still works though)

void resultMatrixToFile(const vector<vector<int>>& matrix);
void matrixToFile(const vector<vector<int>> &A, const vector<vector<int>> &B);
void printMatrix(const vector<vector<int>>& matrix);
void displayMenu();
void testCase();
void divConqMultiplyRecursion(const vector<vector<int>>& A, const vector<vector<int>>& B, vector<vector<int>>& C,
    int rowA, int colA, int rowB, int colB, int rowC, int colC, int size);
vector<vector<int>> generateRandomMatrix(int n, int minVal = -100, int maxVal = 100);
vector<vector<int>> readMatrixFromFile(const string& filename);
vector<vector<int>> addMatrix(const vector<vector<int>> &A, const vector<vector<int>> &B);
vector<vector<int>> subMatrix(const vector<vector<int>> &A, const vector<vector<int>> &B);
vector<vector<int>> multiplyMatrices(const vector<vector<int>>& A, const vector<vector<int>>& B);
vector<vector<int>> divConqMultiply(const vector<vector<int>>& A, const vector<vector<int>>& B);
void strassenMultiplyRecursive(const vector<vector<int>>& A, const vector<vector<int>>& B, 
    vector<vector<int>>& C,
    int rowA, int colA, int rowB, int colB, int rowC, int colC,
    int size, vector<vector<int>>& workspace);
vector<vector<int>> strassenMultiply(const vector<vector<int>>& A, const vector<vector<int>>& B);


// Function to generate a random n x n matrix
vector<vector<int>> generateRandomMatrix(int n, int minVal, int maxVal) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dist(minVal, maxVal);

    vector<vector<int>> matrix(n, vector<int>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matrix[i][j] = dist(gen);
        }
    }
    return matrix;
}

// Function to read a matrix from a file
vector<vector<int>> readMatrixFromFile(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(EXIT_FAILURE);
    }

    int n;
    file >> n;
    vector<vector<int>> matrix(n, vector<int>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!(file >> matrix[i][j])) {
                cerr << "Error: Invalid file format!" << endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    return matrix;
}

// Function to write result into file
void resultMatrixToFile(const vector<vector<int>>& matrix) {
    ofstream file("matrixDatC.txt");  // File name set to matrixDatC.txt
    if (!file) {
        cerr << "Error: Cannot write to file matrixDatC.txt" << endl;
        return;
    }

    int n = matrix.size();
    file << n << endl;
    for (const auto& row : matrix) {
        for (int num : row) {
            file << num << " ";
        }
        file << endl;
    }
}

// Function to write result into file
void matrixToFile(const vector<vector<int>> &A, const vector<vector<int>> &B) {
    static int saveCount = 1;
    string filenameA = "MatrixDatA" + to_string(saveCount) + ".txt";
    string filenameB = "MatrixDatB" + to_string(saveCount) + ".txt";
    ofstream fileA(filenameA); 
    ofstream fileB(filenameB);
    if (!fileA || !fileB) {
        cerr << "Error: Cannot write to file(s)" << endl;
        return;
    }

    int n = A.size();
    fileA << n << endl;
    for (const auto& row : A) {
        for (int num : row) {
            fileA << num << " ";
        }
        fileA << endl;
    }

    fileB << n << endl;
    for (const auto& row : A) {
        for (int num : row) {
            fileB << num << " ";
        }
        fileB << endl;
    }
    cout << "Matrices saved to: " << filenameA <<", "<< filenameB << endl;
}

// Function to print a matrix
void printMatrix(const vector<vector<int>>& matrix) {
    for (const auto& row : matrix) {
        for (int num : row) {
            cout << num << " ";
        }
        cout << endl;
    }
}

//helper function, add matrix 
vector<vector<int>> addMatrix(const vector<vector<int>> &A, const vector<vector<int>> &B){
    int n = A.size();
    vector<vector<int>> C(n, vector<int>(n, 0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = A[i][j] + B[i][j];
    return C;
}

//helper function, sub matrix 
vector<vector<int>> subMatrix(const vector<vector<int>> &A, const vector<vector<int>> &B){
    int n = A.size();
    vector<vector<int>> C(n, vector<int>(n, 0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = A[i][j] - B[i][j];
    return C;

}

// Regular matrix multiplication
vector<vector<int>> multiplyMatrices(const vector<vector<int>>& A, const vector<vector<int>>& B) {
    int n = A.size();
    vector<vector<int>> C(n, vector<int>(n, 0));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

//divide and conquer recursion
void divConqMultiplyRecursion(const vector<vector<int>>& A, const vector<vector<int>>& B, vector<vector<int>>& C,
    int rowA, int colA, int rowB, int colB, int rowC, int colC, int size) {
if (size == 1) {  // Base case
C[rowC][colC] += A[rowA][colA] * B[rowB][colB];
return;
}

int newSize = size / 2;

// Recursively compute submatrices without allocating new ones
divConqMultiplyRecursion(A, B, C, rowA, colA, rowB, colB, rowC, colC, newSize);         // C11 += A11 * B11
divConqMultiplyRecursion(A, B, C, rowA, colA + newSize, rowB + newSize, colB, rowC, colC, newSize); // C11 += A12 * B21

divConqMultiplyRecursion(A, B, C, rowA, colA, rowB, colB + newSize, rowC, colC + newSize, newSize); // C12 += A11 * B12
divConqMultiplyRecursion(A, B, C, rowA, colA + newSize, rowB + newSize, colB + newSize, rowC, colC + newSize, newSize); // C12 += A12 * B22

divConqMultiplyRecursion(A, B, C, rowA + newSize, colA, rowB, colB, rowC + newSize, colC, newSize); // C21 += A21 * B11
divConqMultiplyRecursion(A, B, C, rowA + newSize, colA + newSize, rowB + newSize, colB, rowC + newSize, colC, newSize); // C21 += A22 * B21

divConqMultiplyRecursion(A, B, C, rowA + newSize, colA, rowB, colB + newSize, rowC + newSize, colC + newSize, newSize); // C22 += A21 * B12
divConqMultiplyRecursion(A, B, C, rowA + newSize, colA + newSize, rowB + newSize, colB + newSize, rowC + newSize, colC + newSize, newSize); // C22 += A22 * B22
}

// Wrapper function i hate c++ why do i have to optimize EVERYTHING I HATE IT why can i NEVER avoid POINTERS and when i try to avoid them EVERYTHING IS SO HARD!!!!
vector<vector<int>> divConqMultiply(const vector<vector<int>>& A, const vector<vector<int>>& B) {
int n = A.size();
vector<vector<int>> C(n, vector<int>(n, 0));
divConqMultiplyRecursion(A, B, C, 0, 0, 0, 0, 0, 0, n);
return C;
}

//strassen mult recursive
void strassenMultiplyRecursive(const vector<vector<int>>& A, const vector<vector<int>>& B, 
    vector<vector<int>>& C,
    int rowA, int colA, int rowB, int colB, int rowC, int colC,
    int size, vector<vector<int>>& workspace) {
// Base case
if (size <= 64) {
for (int i = 0; i < size; i++) {
for (int j = 0; j < size; j++) {
int sum = 0;
for (int k = 0; k < size; k++) {
sum += A[rowA + i][colA + k] * B[rowB + k][colB + j];
}
C[rowC + i][colC + j] = sum;
}
}
return;
}

int half = size / 2;

// Indices for workspace (M1-M7)
const int M1 = 0, M2 = 1, M3 = 2, M4 = 3, M5 = 4, M6 = 5, M7 = 6;

// Compute M1 = (A11 + A22) × (B11 + B22)
for (int i = 0; i < half; i++) {
for (int j = 0; j < half; j++) {
int a = A[rowA + i][colA + j] + A[rowA + half + i][colA + half + j];
int b = B[rowB + i][colB + j] + B[rowB + half + i][colB + half + j];
workspace[M1][i * half + j] = a * b;
}
}

// Compute M2 = (A21 + A22) × B11
for (int i = 0; i < half; i++) {
for (int j = 0; j < half; j++) {
int a = A[rowA + half + i][colA + j] + A[rowA + half + i][colA + half + j];
workspace[M2][i * half + j] = a * B[rowB + i][colB + j];
}
}

// Compute M3 = A11 × (B12 - B22)
for (int i = 0; i < half; i++) {
for (int j = 0; j < half; j++) {
int b = B[rowB + i][colB + half + j] - B[rowB + half + i][colB + half + j];
workspace[M3][i * half + j] = A[rowA + i][colA + j] * b;
}
}

// Compute M4 = A22 × (B21 - B11)
for (int i = 0; i < half; i++) {
for (int j = 0; j < half; j++) {
int b = B[rowB + half + i][colB + j] - B[rowB + i][colB + j];
workspace[M4][i * half + j] = A[rowA + half + i][colA + half + j] * b;
}
}

// Compute M5 = (A11 + A12) × B22
for (int i = 0; i < half; i++) {
for (int j = 0; j < half; j++) {
int a = A[rowA + i][colA + j] + A[rowA + i][colA + half + j];
workspace[M5][i * half + j] = a * B[rowB + half + i][colB + half + j];
}
}

// Compute M6 = (A21 - A11) × (B11 + B12)
for (int i = 0; i < half; i++) {
for (int j = 0; j < half; j++) {
int a = A[rowA + half + i][colA + j] - A[rowA + i][colA + j];
int b = B[rowB + i][colB + j] + B[rowB + i][colB + half + j];
workspace[M6][i * half + j] = a * b;
}
}

// Compute M7 = (A12 - A22) × (B21 + B22)
for (int i = 0; i < half; i++) {
for (int j = 0; j < half; j++) {
int a = A[rowA + i][colA + half + j] - A[rowA + half + i][colA + half + j];
int b = B[rowB + half + i][colB + j] + B[rowB + half + i][colB + half + j];
workspace[M7][i * half + j] = a * b;
}
}

// Combine results into C
for (int i = 0; i < half; i++) {
for (int j = 0; j < half; j++) {
// C11 = M1 + M4 - M5 + M7
C[rowC + i][colC + j] = workspace[M1][i * half + j] + workspace[M4][i * half + j] 
         - workspace[M5][i * half + j] + workspace[M7][i * half + j];

// C12 = M3 + M5
C[rowC + i][colC + half + j] = workspace[M3][i * half + j] + workspace[M5][i * half + j];

// C21 = M2 + M4
C[rowC + half + i][colC + j] = workspace[M2][i * half + j] + workspace[M4][i * half + j];

// C22 = M1 - M2 + M3 + M6
C[rowC + half + i][colC + half + j] = workspace[M1][i * half + j] - workspace[M2][i * half + j] 
                     + workspace[M3][i * half + j] + workspace[M6][i * half + j];
}
}
}

//wrapper wrapper wrapper strassen who even are you man
vector<vector<int>> strassenMultiply(const vector<vector<int>>& A, const vector<vector<int>>& B) {
    int n = A.size();
    vector<vector<int>> C(n, vector<int>(n, 0));
    
    // Single pre-allocation of all needed memory
    vector<vector<int>> workspace(7, vector<int>(n * n / 4)); // Enough space for all M matrices
    
    // Start recursion
    strassenMultiplyRecursive(A, B, C, 0, 0, 0, 0, 0, 0, n, workspace);
    return C;
}

// Testing function 
void testCase() {
    int n;
    cout << "Enter matrix size (must be power of 2): ";
    cin >> n;

    cout << "Running 10 trials for " << n << "x" << n << " matrices...\n";

    double total_standard = 0.0, total_strassen = 0.0, total_div_conq = 0.0;

    for (int trial = 1; trial <= 10; trial++) {
        auto A = generateRandomMatrix(n);
        auto B = generateRandomMatrix(n);

        // Standard Multiplication
        auto start = high_resolution_clock::now();
        auto C1 = multiplyMatrices(A, B);
        auto end = high_resolution_clock::now();
        double standard_time = duration_cast<milliseconds>(end - start).count();
        total_standard += standard_time;

        // Divide & Conquer Multiplication
        start = high_resolution_clock::now();
        auto C2 = divConqMultiply(A, B);
        end = high_resolution_clock::now();
        double div_conq_time = duration_cast<milliseconds>(end - start).count();
        total_div_conq += div_conq_time;

        // Strassen Multiplication
        start = high_resolution_clock::now();
        auto C3 = strassenMultiply(A, B);
        end = high_resolution_clock::now();
        double strassen_time = duration_cast<milliseconds>(end - start).count();

        total_strassen += strassen_time;

        //print runtime of each trial
        cout << fixed << setprecision(3);
        
        cout << "Trial " << trial << " - Standard: ";
        if (standard_time < 1000)
            cout << standard_time << "ms, ";
        else
            cout << fixed << setprecision(3) << (standard_time / 1000.0) << "s, ";
        
        cout << "Divide & Conquer: ";
        if (div_conq_time < 1000)
            cout << div_conq_time << "ms, ";
        else
            cout << fixed << setprecision(3) << (div_conq_time / 1000.0) << "s, ";
        
        cout << "Strassen: ";
        if (strassen_time < 1000)
            cout << strassen_time << "ms\n";
        else
            cout << fixed << setprecision(3) << (strassen_time / 1000.0) << "s"<< endl;        
    }

}

// Function to display menu
void displayMenu() {
    cout << "--------- Matrix Multiplication ---------"<< endl;
    cout << "[1]. Naive Matrix Multiplication"<<endl;
    cout << "[2]. Strassen's Matrix Multiplication"<<endl;
    cout << "[3]. Divide and Conquer Matrix Multipliaction"<<endl;
    cout << "[4]. Generate Random Matrices"<<endl;
    cout << "[5]. Load Matrices from Files"<<endl;
    cout << "Filename for matrix A must be 'matrixDatA.txt', B must be 'matrixDatB.txt" << endl;
    cout << "[6]. Save loaded matrices to file, filename 'MatrixAX.txt', 'MatrixDatBX.txt'" << endl;
    cout << "[7]. Run test cases for each algorithm using random data" << endl;
    cout << "[8]. Exit"<<endl;
    cout << "Enter your choice: ";
}

int main() {
    int choice;
    int n;
    vector<vector<int>> A, B, C;
    high_resolution_clock::time_point start, stop;

    while (true) {
        displayMenu();
        cin >> choice;

        switch (choice) {
            case 1: // Naive Multiplication
            if (A.empty() || B.empty()) {
                cout << "Error: Matrices are not initialized!"<< endl;
                break;
            }
            start = high_resolution_clock::now();

            C = multiplyMatrices(A, B);
            cout << "Resultant Matrix (A * B):"<< endl;
            printMatrix(C);

            stop = high_resolution_clock::now();
            cout << "Runtime: " << duration_cast<milliseconds>(stop - start).count() << "ms" << endl;
            resultMatrixToFile(C);  //writes result into file
            cout << "Resultant has been written to matrixDatC.txt." << endl;
            cout << "Matrix Multiplication method"<< endl;
            cout << endl;
            cout << "Press any key to continue." << endl;
            cin.ignore(); 
            cin.get();    
            break;

            case 2: // Strassen’s Multiplication
            if (A.empty() || B.empty()) {
                cout << "Error: Matrices are not initialized!"<< endl;
                break;
                }
                start = high_resolution_clock::now();

            C = strassenMultiply(A, B);
            cout << "Resultant Matrix (A * B):"<< endl;
            printMatrix(C);

            stop = high_resolution_clock::now();
            cout << "Runtime: " << duration_cast<milliseconds>(stop - start).count() << "ms" << endl;
            resultMatrixToFile(C);  //writes result into file
            cout << "Resultant has been written to matrixDatC.txt." << endl;
            cout << "Strassen's Multiplication method"<< endl;
            cout << endl;

            cout << "Press any key to continue." << endl;
            cin.ignore(); 
            cin.get();   
            break;

            case 3: // Divide and Conquer Matrix Multiplication
            if (A.empty() || B.empty()) {
                cout << "Error: Matrices are not initialized!"<< endl;
                break;
            }
            start = high_resolution_clock::now();

            C = divConqMultiply(A, B);
            cout << "Resultant Matrix (A * B):"<< endl;
            printMatrix(C);

            stop = high_resolution_clock::now();
            cout << "Runtime: " << duration_cast<milliseconds>(stop - start).count() << "ms" << endl;

            resultMatrixToFile(C);  //writes result into file
            cout << "Resultant has been written to matrixDatC.txt." << endl;
            cout << "Divide and Conquer method"<< endl;
            cout << endl;
            cout << "Press any key to continue." << endl;

            cin.ignore(); 
            cin.get();    
            break;

            case 4: // Generate Random Matrices
            cout << "Enter matrix size (power of 2): "<< endl;
            cin >> n;
            A = generateRandomMatrix(n);
            B = generateRandomMatrix(n);

            cout << "Matrix A:"<< endl;
            printMatrix(A);
            cout << endl;

            cout << "Matrix B:"<< endl;
            printMatrix(B);
            cout << endl;

            cout << "Press any key to continue." << endl;
            cout << endl;
            cin.ignore(); 
            cin.get();    
            break;

            case 5: // Load Matrices from Files
            A = readMatrixFromFile("matrixDatA.txt");
            B = readMatrixFromFile("matrixDatB.txt");

            cout << "Matrix A (from file):"<< endl;
            printMatrix(A);
            cout << "Matrix B (from file):"<< endl;
            printMatrix(B);
                
            cout << "Press any key to continue." << endl;
            cout << endl;

            cin.ignore(); 
            cin.get();    
            break;
            case 6:
            matrixToFile(A, B);
            cout << "Press any key to continue." << endl;
            cout << endl;

            cin.ignore(); 
            cin.get();    
            break;
            
            case 7:
            testCase();
            cout << "Press any key to continue." << endl;
            cout << endl;

            cin.ignore(); 
            cin.get();    
            break;

            case 8: // Exit
                cout << "Exiting program..."<< endl;
                cout << endl;

                return 0;

            default:
                cout << "Invalid choice! Please try again."<< endl;
        }
    }

    return 0;
}