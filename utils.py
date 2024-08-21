import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def each_data(data):
    for i in range(data.shape[1]):
        plt.figure()

        # Plotting each feature
        sns.lineplot(x=np.arange(data.shape[0]), y=data[:, i])

        # You can also use plt.plot() if you prefer not using seaborn
        # plt.plot(data[:, i])

        plt.title(str(i))
        plt.xlabel('Time Step')
        plt.ylabel('Value')

        # Display plot
        plt.show()

if __name__ == '__main__':
    data_path = 'data/'
    # Load data
    # train test split
    test_list = ['40_60_', '50_10_', '50_20_', '50_30_', '50_40_', '50_50_', '50_60_', '50_70_', '50_80_',
                 '60_10_', '60_20_', '60_30_', '60_40_', '60_50_', '60_60_', '60_70_', '60_80_']
    train_list = []

    for first in range(0, 90, 10):
        for second in range(0, 90, 10):
            orientation = str(first) + '_' + str(second) + '_'
            if orientation not in test_list:
                # test_list.append(orientation)
                train_list.append(orientation)
    print(len(train_list))

    train_list = ['0_0_']

    for orientation in train_list:
        Fp = np.loadtxt(data_path + orientation + 'Fp2000.npy')
        F = np.loadtxt(data_path + orientation + 'F2000.npy')

        print(Fp.shape, F.shape)
        tlen = int(len(Fp) / 8 / 8)
        print(tlen)
        Fp = Fp.reshape(tlen, 8, 8, 9)
        F = F.reshape(tlen, 8, 8, 9)
        Fp = Fp[:, 0, 0, :]
        F = F[:, 0, 0, :]

        print(Fp.shape, F.shape)
        each_data(Fp)
        each_data(F)

