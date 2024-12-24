#include <cstddef>
#include <unordered_map>
#include <iostream>
#include <utility>
#include <cmath>
#include <stdexcept>
#include <algorithm>

template <typename T>
class SparseMatrix {
private:
    struct pair_hash {
        template <class T1, class T2>
        std::size_t operator()(const std::pair<T1, T2>& pair) const {
            std::size_t h1 = std::hash<T1>{}(pair.first);
            std::size_t h2 = std::hash<T2>{}(pair.second);
            return h1 ^ (h2 << 1);
        }
    };

    std::unordered_map<std::pair<size_t, size_t>, T, pair_hash> data;
    size_t rows, cols;

public:
    SparseMatrix(size_t r, size_t c) : rows(r), cols(c) {}

    void AddValue(size_t row, size_t col, T value) {
        auto key = std::make_pair(row, col);
        if (value != T{}) {
            data[key] = value;
        } else {
            data.erase(key);
        }
    }

    T GetValueByIndex(size_t row, size_t col) const {
        auto it = data.find({row, col});
        return (it != data.end()) ? it->second : T{};
    }

    size_t getRows() const { return rows; }
    size_t getCols() const { return cols; }

    std::pair<size_t, size_t> GetSize() const {
        return {rows, cols};
    }

    void Transpose() {
        std::unordered_map<std::pair<size_t, size_t>, T, pair_hash> newData;
        for (const auto& [key, value] : data) {
            newData[{key.second, key.first}] = value;
        }
        data = std::move(newData);
        std::swap(rows, cols);
    }

    void AddMatrix(const SparseMatrix<T>& second_matrix) {
        if (rows != second_matrix.rows || cols != second_matrix.cols) {
            throw std::invalid_argument("Matrix addition is not possible!");
        }

        for (auto& [key, value] : second_matrix.data) {
            data[key] += value;
            if (data[key] == T{}) {
                data.erase(key);
            }
        }
    }

    SparseMatrix MultiplyMatrix(const SparseMatrix<T>& second_matrix) const {
        if (cols != second_matrix.rows) {
            throw std::invalid_argument("Matrix multiplication is not possible!");
        }

        SparseMatrix<T> result(rows, second_matrix.cols);

        // Преобразуем второй аргумент в удобную для поиска структуру
        std::unordered_map<size_t, std::unordered_map<size_t, T>> col_map;
        for (const auto& [key, value] : second_matrix.data) {
            col_map[key.first][key.second] = value;
        }

        for (const auto& [key1, value1] : data) {
            size_t row1 = key1.first;
            size_t col1 = key1.second;

            auto it = col_map.find(col1);
            if (it != col_map.end()) {
                for (const auto& [col2, value2] : it->second) {
                    result.AddValue(row1, col2, result.GetValueByIndex(row1, col2) + value1 * value2);
                }
            }
        }

        return result;
    }

    void PowerMatrix(int exponent) {
        if (rows != cols) {
            throw std::invalid_argument("Only square matrices can be raised to a power.");
        }

        if (exponent == 0) {
            SparseMatrix<T> identity(rows, cols);
            for (size_t i = 0; i < rows; ++i) {
                identity.AddValue(i, i, T(1));
            }
            *this = identity;
            return;
        }

        SparseMatrix<T> result = *this;
        SparseMatrix<T> base = *this;

        for (int exp = 1; exp < exponent; ++exp) {
            result = result.MultiplyMatrix(base);
        }

        *this = result;
    }

    void PrintMatrix() const {
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                std::cout << GetValueByIndex(i, j) << " ";
            }
            std::cout << "\n";
        }
    }

    void AddScalar(T scalar) {
    for (auto& [key, value] : data) {
        value += scalar;
    }
}

void SubtractScalar(T scalar) {
    for (auto& [key, value] : data) {
        value -= scalar;
    }
}

void MultiplyScalar(T scalar) {
    if (scalar == T{}) {
        data.clear();
        return;
    }
    for (auto& [key, value] : data) {
        value *= scalar;
    }
}

void DivideScalar(T scalar) {
    if (scalar == T{}) {
        throw std::invalid_argument("Division by zero.");
    }
    for (auto& [key, value] : data) {
        value /= scalar;
    }
}

void ElementwisePower(T exponent) {
    for (auto& [key, value] : data) {
        value = std::pow(value, exponent);
    }
}

};
template <typename T>
class SparseVector {
private:
    std::unordered_map<size_t, T> data;
    size_t size;

public:
    SparseVector(size_t n) : size(n) {}

    void AddValue(size_t index, T value) {
        if (value != T{}) {
            data[index] = value;
        } else {
            data.erase(index);
        }
    }

    T GetValueByIndex(size_t index) const {
        auto it = data.find(index);
        return (it != data.end()) ? it->second : T{};
    }

    size_t GetSize() const { return size; }

    void PrintElements() const {
        for (const auto& [index, value] : data) {
            std::cout << "Index: " << index << ", Value: " << value << "\n";
        }
    }

    void AddVector(SparseVector& second_vector) {
        if (size != second_vector.size) {
            throw std::invalid_argument("Vector addition is impossible!");
        }

        for (auto& [index, value] : second_vector.data) {
            data[index] += value;
            if (data[index] == T{}) {
                data.erase(index);
            }
        }
    }

    void MultiplyVector(SparseVector& second_vector) {
        if (size != second_vector.size) {
            throw std::invalid_argument("Vector multiplication is impossible!");
        }

        double result = 0;

        for (const auto& [index1, value1] : data) {
            double value2 = second_vector.GetValueByIndex(index1);
            result += value1 * value2;
        }

    }

    SparseVector<T> MultiplyVectorAndMatrix(const SparseMatrix<T>& matrix) {
    if (size != matrix.getRows()) {
        throw std::invalid_argument("Vector and matrix multiplication is not possible!");
    }

    std::unordered_map<size_t, T> result_data;

    // Умножаем вектор на матрицу
    for (typename std::unordered_map<size_t, T>::const_iterator it = data.begin(); it != data.end(); ++it) {
        for (size_t col = 0; col < matrix.getCols(); ++col) {
            T matrix_value = matrix.GetValueByIndex(it->first, col);
            if (matrix_value != T{}) {
                result_data[col] += it->second * matrix_value;
            }
        }
    }

    // Создаем новый результат
    SparseVector<T> result(matrix.getCols());
    for (typename std::unordered_map<size_t, T>::iterator it = result_data.begin(); it != result_data.end(); ++it) {
        if (it->second != T{}) {
            result.AddValue(it->first, it->second);
        }
    }

    return result;  // Возвращаем результат
}


    void Transpose() {
        std::unordered_map<size_t, T> newData;
        for (auto& [index, value] : data) {
            newData[size - 1 - index] = value;
        }
        data = std::move(newData);
    }

    void AddScalar(T scalar) {
        for (auto& [index, value] : data) {
            value += scalar;
        }
    }

    void SubtractScalar(T scalar) {
        for (auto& [index, value] : data) {
            value -= scalar;
        }
    }

    void MultiplyVectorScalar(T scalar) {
        for (auto& [index, value] : data) {
            value *= scalar;
        }
    }

    void DivideScalar(T scalar) {
        if (scalar == T{}) {
            throw std::invalid_argument("Division by zero.");
        }
        for (auto& [index, value] : data) {
            value /= scalar;
        }
    }

    void ElementwisePower(T exponent) {
        for (auto& [index, value] : data) {
            value = std::pow(value, exponent);
        }
    }
};