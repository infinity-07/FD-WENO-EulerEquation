# 编译器
CC = mpicxx

# 编译选项
CFLAGS = -Wall -O4 -g
DEPFLAGS = -MMD -MP

# 目录路径
SRC_DIR = ./src
INC_DIR = ./include
BUILD_DIR = ./build
BIN_DIR = ./bin

# 源文件
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))
DEPS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.d,$(SRCS))

# 可执行文件
TARGET = $(BIN_DIR)/program

# 确保目录存在
$(shell mkdir -p $(BUILD_DIR) $(BIN_DIR))

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) $(DEPFLAGS) -c $< -o $@ -I$(INC_DIR)

-include $(DEPS)

clean:
	rm -f $(OBJS) $(DEPS) $(TARGET)  