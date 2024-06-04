import matplotlib.pyplot as plt
import pandas as pd
# import matplotlib
# matplotlib.use('Agg')
# open data.txt and plot the real values
with open('data.csv') as f:
    # data = f.read()
    # data = data.split('\n')
    data = pd.read_csv('data.csv') 
    i = []
    real = []
    imag = []
    # get the first column of data and plot it
    # load the csv and plot column 1 on the x axis, column 2 on the y axis
    for index, row in data.iterrows():
        i.append(row['index'])
        real.append(row['real'])
        imag.append(row['imag'])

    
    plt.plot(i, real)
    # add a different colored line for the imaginary part
    plt.plot(i, imag, 'r')

    # because it cannot show the plot, save it to a file
    plt.savefig('plot.png')

