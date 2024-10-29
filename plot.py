import matplotlib.pyplot as plt

# Load the data
data = []
with open("Q.dat", "r") as file:
    for line in file:
        x, y = map(float, line.split())
        data.append((x, y))

# Separate data into x and y lists
x_data, y_data = zip(*data)

# Plot
plt.figure()
plt.plot(x_data, y_data, marker='o', linestyle='-', color='b')
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.title("XY Plot from Q.dat")
plt.grid(True)
plt.show()
