import matplotlib.pyplot as plt
  
x1 = [0.02, 0.03, 0.04, 0.05, 0.06]
y1 = [0.0027, 0.0082, 0.0242, 0.0462, 0.0792]
plt.plot(x1, y1, label = "RM(5,2), L=1", marker='o', markerfacecolor='blue', markersize=5)
  
# y2 = [0.0021, 0.007, 0.0218, 0.0455, 0.0764]
# plt.plot(x1, y2, label = "RM(5,2), L=2", marker='+', markerfacecolor='blue', markersize=10)
  
y3 = [0.0023, 0.0082, 0.0212, 0.045, 0.0759]
plt.plot(x1, y3, label = "RM(5,2), L=4", marker='x', markerfacecolor='blue', markersize=10)

y4 = [0.0017, 0.0076, 0.0214, 0.0452, 0.0744]
plt.plot(x1, y4, label = "RM(5,2), L=16", marker='*', markerfacecolor='blue', markersize=10)

x5 = [0.04, 0.05, 0.06, 0.07]
y5 = [0.0327, 0.0911, 0.1815, 0.3087]
plt.plot(x5, y5, label = "RM(7,3), L=1", marker="s", markerfacecolor='blue', markersize=10)

y6 = [0.0036, 0.0159, 0.0489, 0.1209]
plt.plot(x5, y6, label = "RM(7,3), L=4", marker="p", markerfacecolor='blue', markersize=10)

y7 = [0.002, 0.0105, 0.029, 0.0781]
plt.plot(x5, y7, label = "RM(7,3), L=16", marker="8", markerfacecolor='blue', markersize=10)

plt.xlabel('p')
plt.xlim(0.02, 0.09)
plt.ylabel('FER')
plt.ylim(1e-3, 5e-1)
plt.yscale('log')
# giving a title to my graph
plt.title('SCL BSC')
  
# show a legend on the plot
plt.legend()

plt.grid(color='green', linestyle='--', linewidth=0.5)
  
plt.savefig('BSC.png', bbox_inches='tight')
# function to show the plot
plt.show()
