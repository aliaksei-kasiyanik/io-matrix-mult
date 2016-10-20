#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

namespace {

    void MultSimple(const int *__restrict a, const int *__restrict b, int *__restrict c, int n) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    c[i * n + j] += a[i * n + k] * b[k * n + j];
                }
            }
        }
    }

    void readBlockA(ifstream &in, int block_i, int block_j, int block_size, int n, int *__restrict a) {
        int startA = 8;

        int block_start_position = startA + n * block_size * block_i + block_size * block_j;
        int shiftToNextRow = n - block_size;
        in.seekg(block_start_position, ios_base::beg);

        for (int i = 0; i < block_size; ++i) {
            if (i != 0) {
                in.seekg(shiftToNextRow, ios_base::cur);
            }
            char buffer[block_size];
            in.read(buffer, block_size);
            for (int j = 0; j < block_size; ++j) {
                a[i * block_size + j] = (int) buffer[j];
            }
        }

    }

    void readBlockB(ifstream &in, int block_i, int block_j, int block_size, int n, int *__restrict b) {
        int startB = 16 + n * n;

        int block_start_position = startB + n * block_size * block_i + block_size * block_j;
        int shiftToNextRow = n - block_size;

        in.seekg(block_start_position, ios_base::beg);
        for (int i = 0; i < block_size; ++i) {
            if (i != 0) {
                in.seekg(shiftToNextRow, ios_base::cur);
            }
            char buffer[block_size];
            in.read(buffer, block_size);
            for (int j = 0; j < block_size; ++j) {
                b[i * block_size + j] = buffer[j];
            }
        }
    }

    void writeBlockC(ofstream &out, int block_i, int block_j, int block_size, int n, int *__restrict c) {
        int startC = 8;

        int block_start_position = startC + n * block_size * block_i + block_size * block_j;
        int shiftToNextRow = n - block_size;
        out.seekp(block_start_position, ios_base::beg);

        for (int i = 0; i < block_size; ++i) {
            if (i != 0) {
                out.seekp(shiftToNextRow, ios_base::cur);
            }
            char buffer[block_size];
            for (int j = 0; j < block_size; ++j) {
                buffer[j] = (uint8_t) (c[i * block_size + j] % 256);
            }
            out.write(buffer, block_size);
        }
    }


    void BlockMult(int BLOCK) {
        ifstream in("/Users/akasiyanik/FPMI/Tolstikov/io-matrix-mult/in.bin", ios::in | ios::binary);
        ofstream out("/Users/akasiyanik/FPMI/Tolstikov/io-matrix-mult/out.bin", ios::out | ios::binary);
        if (in.is_open() && out.is_open()) {
            int n, m;
            in.read((char *) &n, sizeof(n));
            in.read((char *) &m, sizeof(m));
            out.write((char *) &n, sizeof(n));
            out.write((char *) &m, sizeof(m));
            int blockCountInRow = (int) ceil(n / BLOCK);
            for (int i = 0; i < blockCountInRow; ++i) {
                for (int j = 0; j < blockCountInRow; ++j) {
                    int *c = new int[BLOCK * BLOCK];
                    //fill with zeros??
                    for (int k = 0; k < blockCountInRow; ++k) {
                        int *a = new int[BLOCK * BLOCK];
                        readBlockA(in, i, k, BLOCK, n, a);
                        int *b = new int[BLOCK * BLOCK];
                        readBlockB(in, k, j, BLOCK, n, b);

                        MultSimple(a, b, c, BLOCK);
                    }
                    writeBlockC(out, i, j, BLOCK, n, c);
                }
            }
        }
    }


    void generateInFile(int n) {
        int m = n;
        ofstream file("/Users/akasiyanik/FPMI/Tolstikov/io-matrix-mult/in.bin", ios::out | ios::binary);
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
                uint8_t num = (uint8_t) (rand() % 10);
                file.write((char *) &num, sizeof(num));
            }

            file.close();
        }
    }

    void printInFile() {
        ifstream file("/Users/akasiyanik/FPMI/Tolstikov/io-matrix-mult/in.bin", ios::in | ios::binary);
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
        ifstream file("/Users/akasiyanik/FPMI/Tolstikov/io-matrix-mult/out.bin", ios::in | ios::binary);
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
    ifstream file("/Users/akasiyanik/FPMI/Tolstikov/io-matrix-mult/in.bin", ios::in | ios::binary);
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
        MultSimple(a, b, c, n);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                cout << c[i * n + j] % 256 << " ";
            }
            cout << endl;
        }

    }
}


int main(int argc, char *argv[]) {
    generateInFile(4);
    printInFile();
    BlockMult(2);
    printOutFile();
    test();
    return 0;
}

