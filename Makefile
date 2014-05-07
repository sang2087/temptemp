CC = g++
SRCS = main.cpp GA.cpp BD.cpp
OBJS = $(SRCS:%.cpp=%.o)
TARGET = TSP

.SUFFIXES: .cpp .o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $@ $(OBJS)

clean:
	rm -f $(OBJS) $(TARGET) $(DEPEND_FILE)

run: $(TARGET)
	./$(TARGET)
#
# ifneq ($(MAKECMDGOALS), clean)
# ifneq ($(MAKECMDGOALS), depend)
# ifneq ($(SRCS),)
# -include $(DEPEND_FILE)
# endif
# endif
# endif
