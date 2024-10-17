.SUFFIXES:
FC=ifort
FCFLAGS = -O2 -g
SRC = detstr.f90

TARGET = detstr

all: $(TARGET)

$(TARGET): $(SRC)
	$(FC) $(FCFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET)