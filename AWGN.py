import matplotlib.pyplot as plt
  
x1 = [1, 2, 3, 4, 5]
y1 = [0.84, 0.56, 0.23, 0.05, 0.0043]
plt.plot(x1, y1, label = "L=1", marker='o', markerfacecolor='blue', markersize=5)
  
x2 = [1, 2, 3, 4, 5]
y2 = [0.7, 0.34, 0.088, 0.0074, 0.0001]
plt.plot(x2, y2, label = "L=2", marker='+', markerfacecolor='blue', markersize=10)
  
x3 = [1, 2, 3, 3.5, 4]
y3 = [0.57, 0.204, 0.0286, 0.0068, 0.0009]
plt.plot(x3, y3, label = "L=4", marker='x', markerfacecolor='blue', markersize=10)

x4 = [1, 2, 3, 3.5, 4]
y4 = [0.355, 0.07715, 0.00355, 0.0005, 5e-5]
plt.plot(x4, y4, label = "L=16", marker='*', markerfacecolor='blue', markersize=10)

x5 = [1, 2, 3, 3.5]
y5 = [0.206, 0.0233, 0.00065, 8.3e-5]
plt.plot(x5, y5, label = "L=64", marker="h", markerfacecolor='blue', markersize=10)

plt.xlabel('Eb/N0 (dB)')
plt.xlim(0, 7)
plt.ylabel('FER')
plt.ylim(1e-5, 1e-0)
plt.yscale('log')
# giving a title to my graph
plt.title('RM(8,3) SCL AWGN')
  
# show a legend on the plot
plt.legend()

plt.grid(color='green', linestyle='--', linewidth=0.5)
  
plt.savefig('AWGN_m8_r3.png', bbox_inches='tight')
# function to show the plot
plt.show()
