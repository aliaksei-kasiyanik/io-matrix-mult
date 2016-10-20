#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

namespace {

//    const string IN = "/Users/akasiyanik/FPMI/Tolstikov/io-matrix-mult/in.bin";
//    const string OUT = "/Users/akasiyanik/FPMI/Tolstikov/io-matrix-mult/out.bin";
    const string IN = "input.bin";
    const string OUT = "output.bin";
    const int BLOCK = 100;

    void MultSimple(const int *__restrict a, const int *__restrict b, int *__restrict c, int n1, int m1, int n2, int m2) {
        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < m2; ++j) {
                for (int k = 0; k < m1; ++k) {
                    c[i * m2 + j] += a[i * m1 + k] * b[k * m2 + j];
                }
            }
        }
    }

    void readBlockA(ifstream &in, int block_i, int block_j, int block_size, int block_n, int block_m, int n,
                    int *__restrict a) {
        int startA = 8;
        int shiftToNextRow = n - block_m;
        int block_start_position = startA + n * block_size * block_i + block_size * block_j;
        in.seekg(block_start_position, ios_base::beg);

        for (int i = 0; i < block_n; ++i) {
            if (i != 0) {
                in.seekg(shiftToNextRow, ios_base::cur);
            }
            char buffer[block_m];
            in.read(buffer, block_m);
            for (int j = 0; j < block_m; ++j) {
                a[i * block_m + j] = (int) buffer[j];
            }
        }

    }

    void readBlockB(ifstream &in, int block_i, int block_j, int block_size, int block_n, int block_m, int n,
                    int *__restrict b) {
        int startB = 16 + n * n;

        int shiftToNextRow = n - block_m;
        int block_start_position = startB + n * block_size * block_i + block_size * block_j;

        in.seekg(block_start_position, ios_base::beg);
        for (int i = 0; i < block_n; ++i) {
            if (i != 0) {
                in.seekg(shiftToNextRow, ios_base::cur);
            }
            char buffer[block_m];
            in.read(buffer, block_m);
            for (int j = 0; j < block_m; ++j) {
                b[i * block_m + j] = buffer[j];
            }
        }
    }

    void writeBlockC(ofstream &out, int block_i, int block_j, int block_size, int block_n, int block_m, int n,
                     int *__restrict c) {
        int startC = 8;

        int shiftToNextRow = n - block_m;

        int block_start_position = startC + n * block_size * block_i + block_size * block_j;
        out.seekp(block_start_position, ios_base::beg);

        for (int i = 0; i < block_n; ++i) {
            if (i != 0) {
                out.seekp(shiftToNextRow, ios_base::cur);
            }
            char buffer[block_m];
            for (int j = 0; j < block_m; ++j) {
                buffer[j] = (uint8_t) (c[i * block_m + j] % 256);
            }
            out.write(buffer, block_m);
        }
    }

    void calculateSizeOfBlock(int n, int BLOCK, int i_block, int j_block, int &n_block, int &m_block) {
        if ((i_block + 1) * BLOCK > n) {
            n_block = n % BLOCK;
        } else {
            n_block = BLOCK;
        }
        if ((j_block + 1) * BLOCK > n) {
            m_block = n % BLOCK;
        } else {
            m_block = BLOCK;
        }
    }

    void BlockMult(int BLOCK) {
        ifstream in(IN, ios::in | ios::binary);
        ofstream out(OUT, ios::out | ios::binary);
        if (in.is_open() && out.is_open()) {
            int n, m;

            in.read((char *) &n, sizeof(n));
            in.read((char *) &m, sizeof(m));
            out.write((char *) &n, sizeof(n));
            out.write((char *) &m, sizeof(m));

            int blockCountInRow = (int) ceil((float) n / BLOCK);
            for (int i = 0; i < blockCountInRow; ++i) {
                for (int j = 0; j < blockCountInRow; ++j) {
                    int c_n, c_m;
                    calculateSizeOfBlock(n, BLOCK, i, j, c_n, c_m);
                    int *c = new int[c_n * c_m];
                    for (int i = 0; i < c_n * c_m; ++i) {
                        c[i] = 0;
                    }
                    for (int k = 0; k < blockCountInRow; ++k) {
                        int a_n, a_m;
                        calculateSizeOfBlock(n, BLOCK, i, k, a_n, a_m);
                        int *a = new int[a_n * a_m];
                        readBlockA(in, i, k, BLOCK, a_n, a_m, n, a);

                        int b_n, b_m;
                        calculateSizeOfBlock(n, BLOCK, k, j, b_n, b_m);
                        int *b = new int[b_n * b_m];
                        readBlockB(in, k, j, BLOCK, b_n, b_m, n, b);

                        MultSimple(a, b, c, a_n, a_m, b_n, b_m);
                        delete[] a;
                        delete[] b;
                    }
                    writeBlockC(out, i, j, BLOCK, c_n, c_m, n, c);
                    delete[] c;
                }
            }
        }
    }


    void generateInFile(int n) {
        int m = n;
        ofstream file(IN, ios::out | ios::binary);
        if (file.is_open()) {
            file.write((char *) &n, sizeof(n));
            file.write((char *) &m, sizeof(m));
            for (int i = 0; i < n * m; i++) {
                uint8_t num = (uint8_t) (rand() % 256);
                file.write((char *) &num, sizeof(num));
            }

            file.write((char *) &n, sizeof(n));
            file.write((char *) &m, sizeof(m));
            for (int i = 0; i < n * m; i++) {
                uint8_t num = (uint8_t) (rand() % 256);
                file.write((char *) &num, sizeof(num));
            }

            file.close();
        }
    }

    void printInFile() {
        ifstream file(IN, ios::in | ios::binary);
        if (file.is_open()) {
            int n, m;

            cout << "MATRIX A" << endl;
            file.read((char *) &n, sizeof(n));
            file.read((char *) &m, sizeof(m));
            uint8_t *a = new uint8_t[n * n];
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    uint8_t num;
                    file.read((char *) &num, sizeof(num));
                    a[i * n + j] = num;
                    cout << unsigned(a[i * n + j]) << " ";
                }
                cout << endl;
            }

            cout << "MATRIX B" << endl;
            file.read((char *) &n, sizeof(n));
            file.read((char *) &m, sizeof(m));
            uint8_t *b = new uint8_t[n * n];
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    uint8_t num;
                    file.read((char *) &num, sizeof(num));
                    b[i * n + j] = num;
                    cout << unsigned(b[i * n + j]) << " ";
                }
                cout << endl;
            }
            file.close();
        }
    }

    void printOutFile() {
        ifstream file(OUT, ios::in | ios::binary);
        if (file.is_open()) {
            int n, m;

            cout << "MATRIX C" << endl;
            file.read((char *) &n, sizeof(n));
            file.read((char *) &m, sizeof(m));
            uint8_t *c = new uint8_t[n * n];
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    uint8_t num;
                    file.read((char *) &num, sizeof(num));
                    c[i * n + j] = num;
                    cout << unsigned(c[i * n + j]) << " ";
                }
                cout << endl;
            }
            file.close();
        }
    }
}

void test() {
    ifstream file(IN, ios::in | ios::binary);
    if (file.is_open()) {
        int n, m;

        cout << "TEST" << endl;
        file.read((char *) &n, sizeof(n));
        file.read((char *) &m, sizeof(m));
        int *a = new int[n * n];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                uint8_t num;
                file.read((char *) &num, sizeof(num));
                a[i * n + j] = num;
            }
        }

        file.read((char *) &n, sizeof(n));
        file.read((char *) &m, sizeof(m));
        int *b = new int[n * n];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                uint8_t num;
                file.read((char *) &num, sizeof(num));
                b[i * n + j] = num;
            }
        }
        int *c = new int[n * n];
        MultSimple(a, b, c, n, n, n, n);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                cout << c[i * n + j] % 256 << " ";
            }
            cout << endl;
        }

    }
}


int main(int argc, char *argv[]) {
    generateInFile(1000);
//    printInFile();
    BlockMult(BLOCK);
//    printOutFile();
//    test();
    return 0;
}

