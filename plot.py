import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
# open data.txt and plot the real values
with open('data.txt') as f:
    data = f.read()
    data = data.split('\n')
    x = []
    y = []
    for row in data: 
        r = row.split(' ')
        if len(r) >=  2:
            x.append(r[0])
            y.append(r[1])
    plt.plot(x, y)
    # remove the x and y axes
    m = max(y)
    plt.axis((0, 12, 0, m))
    # set y limit to 0, 1
    # because it cannot show the plot, save it to a file
    plt.savefig('plot.png')

