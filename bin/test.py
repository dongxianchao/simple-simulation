import numpy as np
import matplotlib.pyplot as plt

def cubic_function(x, a, b, c, d):
    return a + b * x + c * np.power(x, 2) + d * np.power(x, 3)

filename = "vector.txt"

# 初始化一个空列表来存储系数
coefficients = []

with open(filename, 'r') as file:
    for line in file:
        data = line.split()  # 逐行读取并按空格分割
        float_data = [float(item) for item in data]
        
        # 使用每4个元素来形成一个子列表
        for i in range(0, len(float_data), 4):
            if i+3 < len(float_data): # 确保有四个元素
                coefficients.append(float_data[i:i+4])

plt.figure(figsize=(10,10))

interval_length = 5.0100200400801306e-04

for i in range(len(coefficients)):
    x = np.linspace(i * interval_length, (i + 1) * interval_length, 100)
    y = cubic_function(x, coefficients[i][0], coefficients[i][1], coefficients[i][2], coefficients[i][3])
    plt.plot(x, y)

plt.ylim([-30, 2])
plt.title("Cubic Spline Interpolation")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)
plt.show()
