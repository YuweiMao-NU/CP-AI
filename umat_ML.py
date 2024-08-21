import numpy as np

def UMAT_ML(model, Fp_list, F_list):
    # print('feature shape: ', Fp_list.shape, F_list.shape)
    # print(F_list.shape)
    input_Fp = Fp_list
    # i1, i2, i3 = input_Fp[:, 0:1], input_Fp[:, 4:5], input_Fp[:, 8:]
    # input_Fp = np.concatenate((i1, i2, i3), axis=1)

    exogenous_input = F_list
    # d1, d2, d3 = exogenous_input[:, 0:1], exogenous_input[:, 4:5], exogenous_input[:, 8:]
    # exogenous_input = np.concatenate((d1, d2, d3), axis=1)
    features = np.concatenate((input_Fp.reshape(1, -1), exogenous_input.reshape(1, -1)), axis=1)
    # print(features.shape)
    output = model.predict(features)

    def_Fp_tau = output[:, 0] * output[:, 4] * output[:, 8] - output[:, 0] * output[:, 5] * output[:,
                                                                                            7] - output[:,
                                                                                                 3] * output[:,
                                                                                                      1] * output[
                                                                                                           :,
                                                                                                           8] + \
                 output[:, 3] * output[:, 2] * output[:, 7] + output[:, 6] * output[:, 1] * output[:,
                                                                                            5] - output[:,
                                                                                                 6] * output[:,
                                                                                                      2] * output[
                                                                                                           :, 4]
    def_Fp_tau = np.maximum(def_Fp_tau, 0.0001)
    d = def_Fp_tau ** (1.0 / 3.0)
    output = output / d

    Fp_tau = output
    Fp_tau = Fp_tau.reshape(3, 3)
    Fp_tau = Fp_tau.tolist()

    return Fp_tau

