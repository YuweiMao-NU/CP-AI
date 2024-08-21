import numpy as np
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_percentage_error
from sklearn.metrics import mean_squared_error
import pickle

def predict_model():
    data_path = ''
    total_loss = 0
    output_list = []

    input_data_list = []
    target_data_list = []

    Fp = np.loadtxt(data_path + 'Fp1500.npy')
    F = np.loadtxt(data_path + 'F1500.npy')
    print(Fp.shape, F.shape)
    tlen = int(len(Fp) / 8 / 8)
    print(tlen)
    Fp = Fp.reshape(tlen, 8, 8, 9)
    F = F.reshape(tlen, 8, 8, 9)
    for elem in range(8):
        for gn in range(8):
            input_data = Fp[:, elem, gn, :]

            for i in range(init_size, 152):
                target_data = Fp[i:i + output_size, elem, gn, :]

                input_Fp = input_data[i - input_size:i]
                # i1, i2, i3 = input_Fp[:, 0:1], input_Fp[:, 4:5], input_Fp[:, 8:]
                # input_Fp = np.concatenate((i1, i2, i3), axis=1)

                exogenous_input = F[i-input_size:i+1, elem, gn, :]
                # d1, d2, d3 = exogenous_input[:, 0:1], exogenous_input[:, 4:5], exogenous_input[:, 8:]
                # exogenous_input = np.concatenate((d1, d2, d3), axis=1)

                features = np.concatenate((input_Fp.reshape(1, -1), exogenous_input.reshape(1, -1)), axis=1)
                # print(features.shape)
                output = model.predict(features)
                mse = mean_squared_error(target_data.reshape(1, 9), output)
                # print(i, mse, Fp[i, elem, gn, 4], output[0, 4])
                print(i, mse, Fp[i, elem, gn, :], output[0, :])
                input_data = np.concatenate((input_data, output), axis=0)
                # print(elem, gn, i, output)

    input_data_list = np.array(input_data_list)
    target_data_list = np.array(target_data_list)

    # return output_list.reshape(-1, 9)

def load_data(ori_list):
    # load data
    input_data_list = []
    target_data_list = []

    for orientation in ori_list:

        Fp = np.loadtxt(data_path + orientation + 'Fp.npy')
        F = np.loadtxt(data_path + orientation + 'F.npy')

        # print(Fp.shape, F.shape)

        for i in range(input_size, Fp.shape[0] - input_size - output_size):
            target_data = Fp[i:i + output_size, :]

            input_Fp = Fp[i - input_size:i]
            # i1, i2, i3 = input_Fp[:, 0:1], input_Fp[:, 4:5], input_Fp[:, 8:]
            # input_Fp = np.concatenate((i1, i2, i3), axis=1)

            exogenous_input = F[i - input_size:i + 1, :]
            # d1, d2, d3 = exogenous_input[:, 0:1], exogenous_input[:, 4:5], exogenous_input[:, 8:]
            # exogenous_input = np.concatenate((d1, d2, d3), axis=1)
            # print(input_Fp.shape, exogenous_input.shape)

            features = np.concatenate((input_Fp.reshape(1, -1), exogenous_input.reshape(1, -1)), axis=1)
            # print(features.shape, target_data.shape)
            input_data_list.append(features)

            target_data_list.append(target_data)

    # print(input_data_list[0].shape)
    input_data_list = np.array(input_data_list)
    target_data_list = np.array(target_data_list)

    return input_data_list, target_data_list

if __name__ == '__main__':

    init_size = 150
    input_size = 50
    output_size = 1
    feature_size = 9
    exogenous_input_size = input_size+1

    data_path = 'data/'

    # Load data
    # train test split
    # test_list = ['40_60_', '50_10_', '50_20_','50_30_','50_40_','50_50_','50_60_','50_70_','50_80_',
    #              '60_10_', '60_20_','60_30_','60_40_','60_50_','60_60_','60_70_','60_80_']
    # train_list = []
    #
    # for first in range(0, 90, 10):
    #     for second in range(0, 90, 10):
    #         orientation = str(first) + '_' + str(second) + '_'
    #         if orientation not in test_list:
    #             # test_list.append(orientation)
    #             train_list.append(orientation)

    train_list = ['0_0_', '0_20_', '0_30_', '0_40_', '0_50_', '0_60_', '0_70_',
                  '10_0_', '10_10_', '10_80_', '20_0_', '30_0_', '40_0_', '50_0_', '60_0_', '70_0_',
                  '70_30_', '70_40_', '70_50_', '70_60_', '80_30_', '80_40_', '80_50_', '80_60_']
    # train_list = ''
    # load train data
    input_data_list, target_data_list = load_data(train_list)
    print('train input size:', input_data_list.shape, 'target size: ', target_data_list.shape)

    target_data_list = target_data_list.reshape(input_data_list.shape[0], -1)
    input_data_list = input_data_list.reshape(input_data_list.shape[0], -1)
    print(input_data_list.shape, target_data_list.shape)

    # Train and save the model
    model = Ridge()
    model.fit(input_data_list, target_data_list)

    model = pickle.dump(model, open('Ridge.model', 'wb'))

    #model = pickle.load(open('Ridge.model', 'rb'))
    #predict_model()

