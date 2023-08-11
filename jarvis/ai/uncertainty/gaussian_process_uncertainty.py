"""
Code to predict properties and their uncertainty.

ML model used: lgbm
"""
# Ref: https://doi.org/10.1021/acsomega.1c03752
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from jarvis.ai.pkgs.utils import regr_scores
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from sklearn.gaussian_process.kernels import WhiteKernel
from sklearn.gaussian_process.kernels import RationalQuadratic

# from sklearn.gaussian_process.kernels import ExpSineSquared
# from sklearn.gaussian_process.kernels import ConstantKernel as C
# import joblib
from joblib import dump
from collections import OrderedDict


def GaussianProcesses(
    x,
    y,
    jid,
    cv=2,
    n_jobs=-1,
    n_iter=10,
    random_state=508842607,
    scoring="neg_mean_absolute_error",
    prop="exfoliation_energy",
    write_files=True,
):
    """
    Perform Gaussian Processes and determine prediction intervals.

    NOTE: it can deal with ~50 descriptors AT MAX ==> reduction of
          the CFID descriptors is needed
    """
    # TODO: Make writing file in proper python format

    # STEP-2: Splitting the data
    # ***************************
    # 90-10% split for train test

    X_train, X_test, y_train, y_test, jid_train, jid_test = train_test_split(
        x, y, jid, random_state=1, test_size=0.1
    )
    # print ('lenx len y',len(x[0]),len(y))

    # STEP-3: Use a specific ML modeli: GP in this case
    # *************************************************
    #
    # Model definition
    # ================
    # The Kernel can be chanded depending on the property under examination.
    # In this example we use a composite kernel that includes long range
    # effects (RBF), medium term irregularities (RationalQuadratic), and
    # white noise
    # Periodic trends can be added using ExpSineSquared
    # Note: RBF has one parameter for each descriptor used ==> if using a
    #       different number of descriptors, the number of parameters has
    #       to be changed
    RQ1 = 316**2 * RationalQuadratic(
        alpha=0.067,
        length_scale=43.5,
        length_scale_bounds=(1e0, 1e2),
        alpha_bounds=(1e-3, 1e1),
    )
    RQ2 = 316**2 * RationalQuadratic(
        alpha=0.00473,
        length_scale=0.0796,
        length_scale_bounds=(1e-2, 1e0),
        alpha_bounds=(1e-3, 1e1),
    )
    RBF1 = 316**2 * RBF(
        length_scale=[
            16,
            1.48e03,
            9.92e04,
            0.858,
            1.02,
            2.58,
            2.72e03,
            2.22e03,
            0.839,
            0.5,
            6.15,
            130,
            3.9e03,
            3.1,
            0.675,
            6.04,
            1.18e03,
            2.9,
            4.56e03,
            333,
            0.0167,
            797,
            2.37e03,
            1.81e03,
        ]
    )
    White = WhiteKernel(noise_level=0.7, noise_level_bounds=(0.05, 0.8))
    kernel = RQ1 + RQ2 + RBF1 + White

    mid_model = GaussianProcessRegressor(
        kernel=kernel,
        alpha=0.0001,
        n_restarts_optimizer=3,
        normalize_y=True,
        random_state=int(774057201),
    )

    # Model fitting
    # =============
    scaler = StandardScaler().fit(X_train)
    scaler.transform(X_train)
    scaler.transform(X_test)

    mid_model.fit(X_train, y_train)

    dump(mid_model, "model.joblib")
    print("")
    print("Params:")
    print(mid_model.get_params(deep=True))
    print("Initial: ", kernel)
    print("Optimum: ", mid_model.kernel_)
    print(
        "Log-Marginal-Likelihood: ",
        mid_model.log_marginal_likelihood(mid_model.kernel_.theta),
    )
    print("")
    print("")

    # PREDICTIONS and UQ
    # ==================
    pred, sigma = mid_model.predict(X_test, return_std=True)
    print("Model       mae     rmse")
    reg_sc = regr_scores(y_test, pred)
    info = OrderedDict()
    info["MAE_Mid"] = reg_sc["mae"]
    actual = y_test
    print("Mid:", round(reg_sc["mae"], 3), round(reg_sc["rmse"], 3))

    fout2 = open("Intervals1.dat", "w")
    line0 = "#    Jid      Observed       pred_Lower       pred_Mid    "
    line1 = line0 + "RealErr(Mid)     "
    line = line1 + "pred_Upper    PredError    Pred_inBounds\n"
    fout2.write(line)
    sum = 0.0
    count = 0
    MAE_err = 0.0
    for ii in range(len(y_test)):
        ss = float(sigma[ii])
        pp = float(pred[ii])
        sum = sum + ss
        lower = pp - ss
        upper = pp + ss
        err = ss
        real_err = float(abs(y_test[ii] - pp))
        err_err = abs(real_err - err)
        MAE_err = MAE_err + err_err
        if (pp >= lower) and (pp <= upper):
            in_bounds = "True"
            count = count + 1
        else:
            in_bounds = "False"
        line2_0 = str(ii) + " " + jid[ii] + " " + str(y_test[ii]) + " "
        line2_1 = line2_0 + str(lower) + " " + str(pp) + " " + str(real_err)
        line2 = (
            line2_1
            + " "
            + str(upper)
            + " "
            + str(ss)
            + " "
            + str(in_bounds)
            + "\n"
        )
        fout2.write(line2)
    print("")
    print("Number of test materials= " + str(len(actual)))
    print(
        "Percentage of in-bound results= "
        + str((float(count) / (len(actual))) * 100)
        + "%"
    )
    print(" ")
    MAE_error = float(MAE_err) / (len(actual))
    print("MAE predicted error = " + str(MAE_error))
    info["MAE_Error"] = MAE_error
    return info
