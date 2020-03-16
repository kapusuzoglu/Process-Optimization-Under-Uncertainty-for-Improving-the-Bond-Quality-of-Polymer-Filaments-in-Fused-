import pandas as pd
import numpy as np
import theano.tensor as tt


# np.random.seed(1234)


def RK_k_eval(ao_init, _gamma_, _theta_, _eta_):

    pi = tt.constant(np.pi)  # theano constant pi value

    kk_ = (_gamma_ / (3 * ao_init * _eta_ * tt.sqrt(pi))) * ((pi - _theta_) * tt.cos(_theta_) + tt.sin(_theta_)) * (
            (pi - _theta_ + tt.sin(_theta_) * (tt.cos(_theta_))) **
            (1 / 2)) / (((pi - _theta_) ** 2) * ((tt.sin(_theta_)) ** 2))

    return kk_


def neck_growth_model(b2, int_temp, bond_model_param, t_final, dt):
    # Neck growth calculation using the experiment temperature
    # == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == %
    # automatic step size Runge - Kutta - Fehlbergintegration method
    # (Burden and Faires, 1985)
    # To overcome numerical instabilities when theta = 0, the initial BC. is fixed at a time valueslightly different than
    # zero and the corresponding value of theta is determined from Eq. 15
    # They found that the majority of neck growth and sintering in the ABS fibers occurred when the interphase temperature
    # was above 200°C, which meant(based on heat transfer analysis and confirmed by experiments) that the nozzle temperature
    # had a large effect on sintering while environment temperature had little effect. In which delta is time step.
    # For this case, delta is set equal to 2 * dt.  dt is the time step that is used for the interval loop
    # == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == %

    # Filament dimensions
    w = 0.8e-3  # Layer Thickness(m)

    # initial radius: half of layer width in meters
    ao = w / 2

    # Material properties - ABS
    # Surface tension
    gamma_r = bond_model_param[0]

    # with a temp.dependent of Delta Gamma / Delta T = - 0.00345 N / m / K
    delta_gamma = -bond_model_param[1]

    # print(T.ones(T.ceil(t_final / dt) + 1),type(T.ones(T.ceil(t_final / dt) + 1)))
    # delta_gamma = -b1 * T.ones(T.ceil(t_final / dt) + 1)

    # Temperature dependent viscosity eta_r = 5100 Viscosity at temp 240 celc
    eta_r = bond_model_param[2]

    # b2 = 0.056; # model parameter for temp dependent viscosity
    # b2 = beta # model parameter for temp dependent viscosity

    kelv = 273.15  # Kelvin conversion

    # idx = T.argmin(T.abs(int_temp - (240 + kelv)))

    # if int_temp[0] > 240 + kelv:
    #     delta_gamma[0: idx] = -0.0005

    t_r_c = 240  # reference temperature in C
    t_r = t_r_c + kelv  # in K

    eta_ = eta_r * tt.exp(-b2 * (int_temp[1] - t_r))

    gamma = gamma_r + delta_gamma * (int_temp[1] - t_r)
    theta = tt.sqrt(2 * dt * gamma / (eta_ * ao))  # Eq. 15

    # pi = tt.constant(np.pi)  # theano constant pi value

    for jjj in range(3, int(t_final / dt) - 1, 2):
        delta_t = 2 * dt  # 2 * time step

        # k1 calculation at t_bond(i / 2)
        eta_1 = eta_r * tt.exp(-b2 * (int_temp[jjj - 1] - t_r))
        gamma_1 = gamma_r + delta_gamma * (int_temp[jjj - 1] - t_r)
        theta_1 = theta
        k1 = RK_k_eval(ao, gamma_1, theta_1, eta_1)

        # k2 calculation
        eta_2 = eta_r * tt.exp(-b2 * (int_temp[jjj] - t_r))
        gamma_2 = gamma_r + delta_gamma * (int_temp[jjj] - t_r)
        theta_2 = theta + dt * k1
        k2 = RK_k_eval(ao, gamma_2, theta_2, eta_2)

        # k3 calculation
        eta_3 = eta_2
        gamma_3 = gamma_2
        theta_3 = theta + dt * k2
        k3 = RK_k_eval(ao, gamma_3, theta_3, eta_3)

        # k4 calculation
        eta_4 = eta_r * tt.exp(-b2 * (int_temp[jjj + 1] - t_r))
        gamma_4 = gamma_r + delta_gamma * (int_temp[jjj + 1] - t_r)
        theta_4 = theta + 2 * dt * k3
        k4 = RK_k_eval(ao, gamma_4, theta_4, eta_4)

        # theta
        theta += (1 / 6) * delta_t * (k1 + 2 * k2 + 2 * k3 + k4)

    # theta_final = theta[-1] # final value of theta to get the final bond length (corresponding the Ref Temp)

    # y = a * sin(theta)
    # save neck growth for each interface on a given layer
    bond_length = 2 * ao * tt.sin(theta) * 1e3  # BL in mm (2*neckRadius)

    return bond_length

# == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == = %
def bond_model(beta_, bondmodelparam, totalparts, interface_temp, t_frame, idx_train):
    # INPUT:
    # bondmodelparam: stores material & model parameters(deterministic)
    # [T_N, v_p, hth, x, y]
    # partinfo: part's information, ...
    # [number of layers, filaments / layer,  total interfaces, and total filaments, T_N v_p hth] in a part

    ctr_ = 0  # counter over number of parts
    dt = 0.1  # time step for the neck growth prediction
    predicted_bl_train = tt.vector()

    # Predict BLfor each interface of each part
    for partid in range(totalparts):  # loop over each manufactured part ID

        int_temp = interface_temp[ctr_]

        # store bond lengths(at Ref_temp) at each interface
        lastnonzeroidx = int_temp.argmin(axis=0)  # return index min. value of temp of each column

        # Predict Neck growth at each interface
        for int_id in idx_train[idx_train < int_temp.shape[1]]:  # loop over each interface
            # for int_id in range(numtotalint[ctr_]): # loop over each interface

            # Neck growth calculation using the temperature field coming from the Abaqus model
            t_final = t_frame[lastnonzeroidx[int_id]]  # final time value

            # Model Bond Length
            bondlengths = neck_growth_model(beta_, int_temp[0:lastnonzeroidx[int_id], int_id], bondmodelparam,
                                            t_final, dt)

            # Get the bond length at a time when Temperature = Ref_Temp, save BL in mm
            predicted_bl_train = tt.concatenate([bondlengths])



        idx_train = idx_train[idx_train >= int_temp.shape[1]] - int_temp.shape[1]
        # np.delete(idx_train, idx_train<int_temp.shape[1]) # remove the interfaces that are looped over (previous part)

        ctr_ += 1  # increment counter of number of parts
    # predicted_bl_train = (predicted_bl_train - tt.mean(predicted_bl_train)) / tt.std(predicted_bl_train)
    return predicted_bl_train


def neck_growth_model_2(b2, int_temp, bond_model_param, t_final, dt):
    # Neck growth calculation using the experiment temperature
    # == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == %
    # automatic step size Runge - Kutta - Fehlbergintegration method
    # (Burden and Faires, 1985)
    # To overcome numerical instabilities when theta = 0, the initial BC. is fixed at a time valueslightly different than
    # zero and the corresponding value of theta is determined from Eq. 15
    # They found that the majority of neck growth and sintering in the ABS fibers occurred when the interphase temperature
    # was above 200°C, which meant(based on heat transfer analysis and confirmed by experiments) that the nozzle temperature
    # had a large effect on sintering while environment temperature had little effect. In which delta is time step.
    # For this case, delta is set equal to 2 * dt.  dt is the time step that is used for the interval loop
    # == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == %

    # Filament dimensions
    w = 0.8e-3  # Layer Thickness(m)

    # initial radius: half of layer width in meters
    ao = w / 2

    # Material properties - ABS
    # Surface tension
    gamma_r = bond_model_param[0]

    # with a temp.dependent of Delta Gamma / Delta T = - 0.00345 N / m / K
    delta_gamma = -bond_model_param[1]

    # print(T.ones(T.ceil(t_final / dt) + 1),type(T.ones(T.ceil(t_final / dt) + 1)))
    # delta_gamma = -b1 * T.ones(T.ceil(t_final / dt) + 1)

    # Temperature dependent viscosity eta_r = 5100 Viscosity at temp 240 celc
    eta_r = bond_model_param[2]

    # b2 = 0.056; # model parameter for temp dependent viscosity
    # b2 = beta # model parameter for temp dependent viscosity

    kelv = 273.15  # Kelvin conversion

    # idx = T.argmin(T.abs(int_temp - (240 + kelv)))

    # if int_temp[0] > 240 + kelv:
    #     delta_gamma[0: idx] = -0.0005

    t_r_c = 240  # reference temperature in C
    t_r = t_r_c + kelv  # in K

    eta_ = eta_r * np.exp(-b2 * (int_temp[1] - t_r))

    gamma = gamma_r + delta_gamma * (int_temp[1] - t_r)
    theta = np.sqrt(2 * dt * gamma / (eta_ * ao))  # Eq. 15

    pi = (np.pi)  # theano constant pi value

    for jjj in range(3, int(t_final / dt) - 1, 2):
        delta_t = 2 * dt  # 2 * time step

        # k1 calculation at t_bond(i / 2)
        eta_1 = eta_r * np.exp(-b2 * (int_temp[jjj - 1] - t_r))
        gamma_1 = gamma_r + delta_gamma * (int_temp[jjj - 1] - t_r)
        theta_1 = theta
        k1 = (gamma_1 / (3 * ao * eta_1 * np.sqrt(pi))) * ((pi - theta_1) * np.cos(theta_1) +
                                                           np.sin(theta_1)) * (
                     (pi - theta_1 + np.sin(theta_1) * (np.cos(theta_1))) **
                     (1 / 2)) / (((pi - theta_1) ** 2) * ((np.sin(theta_1)) ** 2))

        # k2 calculation
        eta_2 = eta_r * np.exp(-b2 * (int_temp[jjj] - t_r))
        gamma_2 = gamma_r + delta_gamma * (int_temp[jjj] - t_r)
        theta_2 = theta + dt * k1
        k2 = (gamma_2 / (3 * ao * eta_2 * np.sqrt(pi))) * ((pi - theta_2) * np.cos(theta_2) +
                                                           np.sin(theta_2)) * (
                     (pi - theta_2 + np.sin(theta_2) * (np.cos(theta_2))) **
                     (1 / 2)) / (((pi - theta_2) ** 2) * ((np.sin(theta_2)) ** 2))

        # k3 calculation
        eta_3 = eta_2
        gamma_3 = gamma_2
        theta_3 = theta + dt * k2
        k3 = (gamma_3 / (3 * ao * eta_3 * np.sqrt(pi))) * ((pi - theta_3) * np.cos(theta_3) +
                                                           np.sin(theta_3)) * (
                     (pi - theta_3 + np.sin(theta_3) * (np.cos(theta_3))) **
                     (1 / 2)) / (((pi - theta_3) ** 2) * ((np.sin(theta_3)) ** 2))

        # k4 calculation
        eta_4 = eta_r * np.exp(-b2 * (int_temp[jjj + 1] - t_r))
        gamma_4 = gamma_r + delta_gamma * (int_temp[jjj + 1] - t_r)
        theta_4 = theta + 2 * dt * k3
        k4 = (gamma_4 / (3 * ao * eta_4 * np.sqrt(pi))) * ((pi - theta_4) * np.cos(theta_4) +
                                                           np.sin(theta_4)) * (
                     (pi - theta_4 + np.sin(theta_4) * (np.cos(theta_4))) **
                     (1 / 2)) / (((pi - theta_4) ** 2) * ((np.sin(theta_4)) ** 2))

        # theta
        theta += (1 / 6) * delta_t * (k1 + 2 * k2 + 2 * k3 + k4)

    # theta_final = theta[-1] # final value of theta to get the final bond length (corresponding the Ref Temp)

    # y = a * sin(theta)
    # save neck growth for each interface on a given layer
    bond_length = 2 * ao * np.sin(theta) * 1e3  # BL in mm (2*neckRadius)

    return bond_length

# == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == = %
def bond_model_2(beta_, bondmodelparam, totalparts, interface_temp, t_frame, idx_train):
    # INPUT:
    # bondmodelparam: stores material & model parameters(deterministic)
    # [T_N, v_p, hth, x, y]
    # partinfo: part's information, ...
    # [number of layers, filaments / layer,  total interfaces, and total filaments, T_N v_p hth] in a part

    ctr_ = 0  # counter over number of parts
    dt = 0.1  # time step for the neck growth prediction
    predicted_bl_train = []

    # Predict BLfor each interface of each part
    for partid in range(totalparts):  # loop over each manufactured part ID

        int_temp = interface_temp[ctr_]

        # store bond lengths(at Ref_temp) at each interface
        lastnonzeroidx = int_temp.argmin(axis=0)  # return index min. value of temp of each column

        # Predict Neck growth at each interface
        for int_id in idx_train[idx_train < int_temp.shape[1]]:  # loop over each interface
            # for int_id in range(numtotalint[ctr_]): # loop over each interface

            # Neck growth calculation using the temperature field coming from the Abaqus model
            t_final = t_frame[lastnonzeroidx[int_id]]  # final time value

            # Model Bond Length
            bondlengths = neck_growth_model_2(beta_, int_temp[0:lastnonzeroidx[int_id], int_id], bondmodelparam,
                                            t_final, dt)

            # Get the bond length at a time when Temperature = Ref_Temp, save BL in mm
            # predicted_bl_train = np.concatenate([bondlengths])
            predicted_bl_train.append(bondlengths)


        idx_train = idx_train[idx_train >= int_temp.shape[1]] - int_temp.shape[1]
        # np.delete(idx_train, idx_train<int_temp.shape[1]) # remove the interfaces that are looped over (previous part)

        ctr_ += 1  # increment counter of number of parts
    # predicted_bl_train = (predicted_bl_train - tt.mean(predicted_bl_train)) / tt.std(predicted_bl_train)
    predicted_bl_train = np.array(predicted_bl_train)
    return predicted_bl_train


def measurements_and_training_data(num_parts, ratio_):
    # =================================================================================================================== %
    # _______________________________________________       INPUTS     ___________________________________________________
    # =================================================================================================================== %
    global time_frame
    wth = 0.8e-3  # width of filaments [m]

    # Reference temperature at which the neck growth is predicted
    ref_temp = 130  # in celcius
    kelv = 273.15  # 0 Celcius in K

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

    # =================================================================================================================== %
    # ___________________________________        IMPORT Measured BOND LENGTH DATA     ____________________________________
    # =================================================================================================================== %
    # Directory where measurements are stored
    dir_measurements = "C:/Users/berkc/Box Sync/AM/SecondPaperMultiObjective/ModelCalibration/ManufacturedParts_BL"
    dir_temp_data = "C:/Users/berkc/Box Sync/AM/SecondPaperMultiObjective/ModelCalibration/AbaqusTempData" \
                    "/TempData_Sample"

    # IMPORT Measured BOND LENGTH DATA
    start_part_id = 1  # starting ID part
    # num_parts = 2                                          # number of parts(the last value)
    total_parts = np.arange(start_part_id, num_parts + 1, dtype='int32')
    index_to_remove = [4, 13, 20]  # remove parts that do not have measurements, i.e., 5,14,21
    total_parts = np.delete(total_parts, index_to_remove)  # remove elements

    ctr = 0  # counter of parts
    # num_total_int = 0                              # counter of interfaces
    int_temperature = []  # list stores interface temperature data of all parts' filaments
    measured_bl_row = []  # list stores bond length (BL) data of all parts
    inp = []
    # part_info = np.zeros(total_parts.shape[0], dtype='int')

    for part_id in total_parts:  # loop over each manufactured part ID

        # Read Process parameters & bond length measurements
        bl_file = dir_measurements + "/PPSet" + str(part_id) + ".xlsx"

        df = pd.read_excel(bl_file, header=None)  # read the excel file

        num_layers = df.iloc[2, 0]  # # of layers
        num_interfaces = df.iloc[-2, 1]  # # of interfaces/layer
        num_filaments = num_interfaces + 1  # # of filaments/layer

        # num_layers = int(num_layers)
        # num_interfaces = int(num_interfaces)

        num_interfaces_of_a_part = int(num_layers * num_interfaces)  # num. of interfaces of that part
        num_filaments_of_a_part = int(num_layers * num_filaments)  # num. of filaments of that part

        # # # of total interface in all parts
        # num_total_int = num_total_int + num_interfaces_of_a_part

        # save printer temperature, speed, height input for each part
        t_n = df.iloc[0, 1]  # Printer nozzle temperature(ºC)
        v_p = df.iloc[0, 3]  # printer speed mm / s
        hth = df.iloc[0, 5]  # layer height[mm]
        t_n = float(t_n)
        v_p = float(v_p)

        # convert data types to float32 for theano
        # hth = hth.astype('float32')

        print('T_N, v_p, height:', t_n, v_p, hth, "\n")

        raw_measuredbl = df.iloc[2:-1, 3]  # measured bond lengths between each interface
        raw_measuredbl = raw_measuredbl.astype('float32')

        # reshape the measured bond length array & convert to numpy ndarray
        reshaped_measured_bl = raw_measuredbl.values.reshape(num_interfaces, num_layers, order='F').copy()

        # first column is 1st layer and soon(each row is each interface bond length, BL)
        measured_bl = np.fliplr(reshaped_measured_bl)  # flip matrix left to right

        # store measured BL data of all parts in order reshaped in row
        measured_bl_row.append([measured_bl.reshape(num_interfaces_of_a_part).copy()])

        # == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == = %
        # Abaqus MODEL
        # open temperature data file for each node and save to time(t) and Temperature in Celcius(T_C)

        # x: length(constant 35 mm)
        # y: width
        # z: height
        zcoord = np.zeros(num_filaments_of_a_part, dtype='float32')
        t_c = []

        # Read Abaqus Temperature data
        folder_name = dir_temp_data + str(part_id)
        file_name = "Job_fdmTher_int" + str(part_id) + '_'
        dir_folder_filename = folder_name + "/" + file_name

        # extract time array once outside the loop
        filename_time = dir_folder_filename + str(1) + ".txt"

        data = pd.read_csv(filename_time, sep=";", header=None, index_col=None)
        data = data.apply(pd.to_numeric, errors='coerce')

        # import time frame data only once since it is same for all filaments
        time_frame = data.iloc[1:, 3]

        k = 0
        # loop over total number of filaments & read Temperature Data
        for ii in range(1, num_filaments_of_a_part + 1):
            filename = dir_folder_filename + str(ii) + ".txt"  # file name for each filament
            data = pd.read_csv(filename, sep=";", header=None)  # read file
            data = data.apply(pd.to_numeric, errors='coerce')
            zcoord[k] = data.iloc[1, 2]  # read z-coordinate value
            t_c.append(data.iloc[1:, 4])  # read temperature data
            k += 1

        # Convert built Python list of strings to a Numpy float array.
        t_c = np.array(t_c, dtype='float32')

        # convert Celcius to Kelvin
        temp = t_c.T + kelv

        # Average two adjacent lines' temperatures to obtain interface temperatures for each interface on that layer
        # initialize matrix for Interface Temperatures
        int_temp = np.zeros((temp.shape[0], num_interfaces_of_a_part), dtype='float32')

        # Interpolate interface temperatures
        kk_int = 0  # count  # of interfaces
        for jj in range(1, num_layers + 1):  # loop over layers
            for ii in range(2, num_filaments + 1):  # loop over filaments on jj - th layer

                kk_fil = (jj - 1) * (num_filaments - 1) + ii + jj - 2  # filament ID

                # find the index value when (kk + 1)th filament is started to extruded
                t_idx = np.argmax(temp[:, kk_fil] < t_n + kelv)  # Since argmaxx will stop at the first True

                # find the first index that temp value is lower than Ref_temp
                t_idx_ = np.argmax((temp[:, kk_fil] + temp[:, kk_fil - 1]) / 2 < ref_temp + kelv) + 1

                # average two filament's temp. to obtain interface temp.
                int_temp[0:t_idx_ - t_idx, kk_int] = (temp[t_idx:t_idx_, kk_fil - 1] + temp[t_idx:t_idx_, kk_fil]) / 2
                kk_int += 1

        # store each parts' interface temperature in cell array
        int_temperature.append(int_temp)

        # Prediction inputs are x & y coordinates of vertical bond length locations
        # x, y coordinate of layer 1 & interface 1(vertical)
        # x = wth;
        ycoord = 0.5 * hth  # 0.5*height of a layer in mm
        iki_y = ycoord * 2

        # store inputs for GP(model disrepancy at each interface)
        for jj in range(num_layers):
            for ii in range(num_interfaces):
                # use x & y coordinates of vertical bonds as training data for the GP
                # Inp =[ Temperature, speed, height, x, y ]
                inp.append([t_n, v_p, hth, ii * wth, ycoord + (jj - 1) * iki_y])

        # Store PartInfo array: [number of total interfaces, hth] in a part
        # part_info.append([num_interfaces_of_a_part, hth])
        # part_info[ctr] = num_interfaces_of_a_part

        ctr += 1  # increment counter to keep track of total number of parts analyzed

        print('Numberofpart', "NumLayers", "NumInterfaces")
        print("\t", ctr, "\t", "\t", "\t", num_layers, "\t", "\t", "\t", num_interfaces)

    # Inp: stored inputs for Gaussian process
    # (model disrepancy at each interface):
    #           [T_N, v_p, hth, x, y]
    # Convert built Python lists to a Numpy array.
    inp = np.array(inp, dtype='float32')
    int_temperature = np.array(int_temperature)

    # concatenating different size arrays stored in a list
    measured_bl_row = np.concatenate(measured_bl_row, axis=1)
    measured_bl_row = measured_bl_row.T  # transpose s.t. the number of rows matches Inp

    # Normalize training data


    # inp_normalized = (inp - mean_of_data) / std_of_data  # take mean of each column
    # measured_bl_row = (measured_bl_row - measured_bl_row.mean(axis=0)) / measured_bl_row.std(axis=0)  # take mean of each column
    # inp_normalized = inp
    alldata = np.hstack((inp, measured_bl_row))  # stack 2 numpy arrays column-wise

    # -------------------------------------------------------------------------
    #               Random Permutation of Training Data
    # -------------------------------------------------------------------------
    nl = inp.shape[0]  # size of training data

    # randomly select RatioToBeUsed to be training set for GP model
    num_train = round(ratio_ * nl)
    idx_ = np.random.permutation(nl)
    # idx_ = np.arange(nl)  # do not do random permutation

    # Use the first RatioToBeUsed to train the model
    idx_train = idx_[0:num_train]
    all_data_train = alldata[idx_train, :]

    # mean_of_data = all_data_train.mean(axis=0)
    # std_of_data = all_data_train.std(axis=0)
    # all_data_train = (all_data_train - mean_of_data) / std_of_data  # take mean of each column


    # The (1-RatioToBeUsed) will be used to test the model
    idx_test = idx_[(num_train + 1):]
    all_data_test = alldata[idx_test, :]

    x = all_data_train[:, :-1]  # training data, for all but last column
    yy = all_data_train[:, -1]  # measurements of the training data, last column

    return x, yy, ctr, time_frame, int_temperature, idx_train
