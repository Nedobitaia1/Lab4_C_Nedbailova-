#include <cstddef>
#include <random>
#include <iostream>
#include <unordered_set>
#include <stdexcept>
#include <chrono>
#include <vector>
#include "Sparse.cpp"


std::vector<int> sumVectors(const std::vector<int> &v1, const std::vector<int> &v2) {
  if (v1.size() != v2.size()) {
    std::cout<< "The vectors are not equal, and addition is impossible.\n";
    return {};
  }
  std::vector<int> result(v1.size());
  for (size_t i = 0; i < v1.size(); i++) {
    result[i] = v1[i] + v2[i];
  }
  return result;
}

int dotVector(const std::vector<int> &v1, const std::vector<int> &v2) {
  if (v1.size() != v2.size()) {
    throw std::runtime_error("Vector sizes do not match for dot product.");
  }
  int result = 0;
  for (size_t i = 0; i < v1.size(); ++i) {
    result += v1[i] * v2[i];
  }
  return result;
}

std::vector<std::vector<int>>
sumMatrix(const std::vector<std::vector<int>> &m1,
               const std::vector<std::vector<int>> &m2) {
  if (m1.size() != m2.size() || m1.empty() || m1[0].size() != m2[0].size()) {
    throw std::runtime_error("Matrix sizes do not match for addition.");
  }
  size_t rows = m1.size();
  size_t cols = m1[0].size();
  std::vector<std::vector<int>> result(rows, std::vector<int>(cols, 0));
  for (size_t r = 0; r < rows; ++r) {
    for (size_t c = 0; c < cols; ++c) {
      result[r][c] = m1[r][c] + m2[r][c];
    }
  }
  return result;
}

std::vector<std::vector<int>>
mulMatrix(const std::vector<std::vector<int>> &m1,
               const std::vector<std::vector<int>> &m2) {
  if (m1.empty() || m2.empty() ||  m1[0].size() != m2.size()) {
    throw std::runtime_error(
        "Matrix dimensions are incompatible for multiplication.");
  }
  size_t rows = m1.size();
  size_t cols = m2[0].size();
  size_t inner = m1[0].size();
  std::vector<std::vector<int>> result(rows, std::vector<int>(cols, 0));
  for (size_t r = 0; r < rows; ++r) {
    for (size_t c = 0; c < cols; ++c) {
      for (size_t k = 0; k < inner; ++k) {
        result[r][c] += m1[r][k] * m2[k][c];
      }
    }
  }
  return result;
}

std::vector<int> mulOnVectorMatrix(const std::vector<std::vector<int>> &m,
                                        const std::vector<int> &v) {
  if (m.empty() || m[0].size() != v.size()) {
    throw std::runtime_error(
        "Matrix and vector sizes do not match for multiplication.");
  }
  size_t rows = m.size();
  size_t cols = m[0].size();
  std::vector<int> result(rows, 0);
  for (size_t r = 0; r < rows; ++r) {
    for (size_t c = 0; c < cols; ++c) {
      result[r] += m[r][c] * v[c];
    }
  }
  return result;
}




int main() {
    const size_t matrix_size = 1000;          
    const size_t num_elements = matrix_size * matrix_size; 
    const size_t num_nonzeros_matrix = num_elements / 100; 

    std::mt19937 gen(1);  
    std::uniform_int_distribution<int> dist_val(1, 100); 
    std::uniform_int_distribution<size_t> dist_idx(0, num_elements - 1); 

    SparseMatrix<int> sm1(matrix_size, matrix_size);
    SparseMatrix<int> sm2(matrix_size, matrix_size);

    std::unordered_set<size_t> unique_indices;
    while (unique_indices.size() < num_nonzeros_matrix) {
        unique_indices.insert(dist_idx(gen));
    }

    for (size_t idx : unique_indices) {
        size_t row = idx / matrix_size;
        size_t col = idx % matrix_size;
        sm1.AddValue(row, col, dist_val(gen));
        sm2.AddValue(row, col, dist_val(gen));
    }

    std::vector<std::vector<int>> dm1(matrix_size, std::vector<int>(matrix_size, 0));
    std::vector<std::vector<int>> dm2(matrix_size, std::vector<int>(matrix_size, 0));

    for (size_t r = 0; r < matrix_size; ++r) {
        for (size_t c = 0; c < matrix_size; ++c) {
            dm1[r][c] = sm1.GetValueByIndex(r, c);
            dm2[r][c] = sm2.GetValueByIndex(r, c);
        }
    }


    auto start = std::chrono::high_resolution_clock::now();
    sm1.AddMatrix(sm2);           
    auto end = std::chrono::high_resolution_clock::now();
    auto d1 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<int>> dm_add = sumMatrix(dm1, dm2);
    end = std::chrono::high_resolution_clock::now();
    auto d2 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "Matrix Addition: "<< "sparse: " << d1 << " ms; "<< "vector: " << d2 << " ms\n";

    start = std::chrono::high_resolution_clock::now();
    SparseMatrix<int> sm_mul = sm1.MultiplyMatrix(sm2);  
    end = std::chrono::high_resolution_clock::now();
    d1 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<int>> dm_mul = mulMatrix(dm1, dm2);
    end = std::chrono::high_resolution_clock::now();
    d2 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "Matrix Multiplication: "<< "sparse: " << d1 << " ms; "<< "vector: " << d2 << " ms\n";

    start = std::chrono::high_resolution_clock::now();
    sm1.Transpose();  
    end = std::chrono::high_resolution_clock::now();
    d1 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<int>> dm_transpose(matrix_size, std::vector<int>(matrix_size, 0));
    for (size_t r = 0; r < matrix_size; ++r) {
        for (size_t c = 0; c < matrix_size; ++c) {
            dm_transpose[c][r] = dm1[r][c];  
        }
    }
    end = std::chrono::high_resolution_clock::now();
    d2 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "Matrix Transposition: "<< "sparse: " << d1 << " ms; "<< "vector: " << d2 << " ms\n";


    const size_t vector_size = 1000000;          
    const size_t num_nonzeros_vector = vector_size / 100; 

    SparseVector<int> sv1(vector_size);
    SparseVector<int> sv2(vector_size);

    std::unordered_set<size_t> unique_indices_v;
    while (unique_indices_v.size() < num_nonzeros_vector) {
        unique_indices_v.insert(dist_idx(gen));
    }

    for (size_t idx : unique_indices_v) {
        sv1.AddValue(idx, dist_val(gen));
        sv2.AddValue(idx, dist_val(gen));
    }

    std::vector<int> dv1(vector_size, 0);
    std::vector<int> dv2(vector_size, 0);

    for (size_t idx : unique_indices_v) {
        dv1[idx] = sv1.GetValueByIndex(idx);
        dv2[idx] = sv2.GetValueByIndex(idx);
    }

    start = std::chrono::high_resolution_clock::now();
    SparseVector<int> sv_add = sv1;
    sv_add.AddVector(sv2);
    end = std::chrono::high_resolution_clock::now();
    d1 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    std::vector<int> dv_add = sumVectors(dv1, dv2);
    end = std::chrono::high_resolution_clock::now();
    d2 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "Vector Addition: "<< "sparse: " << d1 << " ms; "<< "vector: " << d2 << " ms\n";

    start = std::chrono::high_resolution_clock::now();
    sv1.MultiplyVector(sv2);
    end = std::chrono::high_resolution_clock::now();
    d1 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    start = std::chrono::high_resolution_clock::now();
    int dv_dot = dotVector(dv1, dv2);
    end = std::chrono::high_resolution_clock::now();
    d2 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "Vector Dot Product: " << "sparse: " << d1 << " ms; "<< "vector: " << d2 << " ms\n";

    
    return 0;
}
