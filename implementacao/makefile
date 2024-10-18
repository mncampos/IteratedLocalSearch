# Nome: Mateus Nunes Campos

CXX := g++
CXXFLAGS := -std=c++20 -Wall -O2
LDFLAGS := -static

TARGET := iterated_local_search.exe

SRC := iterated_local_search.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

clean:
	del $(TARGET)
