TARGET = main.exe

CXX = C:/msys64/mingw64/bin/g++

CXXFLAGS = -std=c++14 -Wall

# Добавим Logger.cpp в список исходных файлов
SRCS = main.cpp Sparse.cpp

# Генерация объектных файлов
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	del $(OBJS) $(TARGET)
