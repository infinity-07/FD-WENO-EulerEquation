# 编译器
CC = mpicxx

# 编译选项
# CFLAGS = -Wall -O4 -g -fopenmp
CFLAGS = -Wall -O4 -g
DEPFLAGS = -MMD -MP

# 源文件
SRCS = main.cpp WENOFV.cpp 
OBJS = $(SRCS:.cpp=.o)
DEPS = $(SRCS:.cpp=.d)

# 生成可执行文件
TARGET = program

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.cpp
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

-include $(DEPS)

clean:
	rm -f *.o *.d gmon.out $(TARGET)