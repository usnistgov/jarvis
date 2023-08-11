"""
Code to predict properties and their uncertainty.

ML model used: lgbm
"""
# Ref: https://doi.org/10.1021/acsomega.1c03752
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import pickle
from jarvis.ai.pkgs.utils import regr_scores
from collections import OrderedDict
from jarvis.ai.pkgs.lgbm.regression import get_lgbm


def quantile_regr_predint(
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
    Perform Quantile regression and determine prediction intervals.

    LOWER_ALPHA = 0.16
    Mid model uses ls as loss function, not quantile, to
        optimize for the mean, not the median
    UPPER_ALPHA = 0.84
    This choice of LOWER_ALPHA, UPPER_ALPHA gives a prediction
    interval ideally equal to 0.68, i.e. 1 standard deviation.
    However, the number of in-bound prediction must be computed
    for the specific fitted models, and that gives the true meaning
    of the uncertainties computed here.
    See:
    https://machinelearningmastery.com/prediction-intervals-for-machine-learning
    https://www.inovex.de/blog/uncertainty-quantification-deep-learning
    """
    # TODO: Make writing file in proper python format

    # STEP-2: Splitting the data
    # ***************************
    # 90-10% split for train test

    X_train, X_test, y_train, y_test, jid_train, jid_test = train_test_split(
        x, y, jid, random_state=1, test_size=0.1
    )
    # print ('lenx len y',len(x[0]),len(y))

    # STEP-3: Use a specific ML model
    # ********************************

    # Set lower and upper quantile
    # StanDev
    LOWER_ALPHA = 0.16
    # MID_ALPHA = 0.50
    UPPER_ALPHA = 0.84

    # LOWER Model
    # ===========
    scaler = StandardScaler().fit(X_train)
    scaler.transform(X_train)
    scaler.transform(X_test)

    objective = "quantile"
    alpha = LOWER_ALPHA
    print("Prima di lgbm for LOWER model")
    lower_model = get_lgbm(
        X_train,
        X_test,
        y_train,
        y_test,
        cv,
        n_jobs,
        scoring,
        n_iter,
        objective,
        alpha,
        random_state,
    )
    print("Dopo lgbm for LOWER model")
    name = str(prop) + str("_lower")
    filename = str("pickle2-") + str(name) + str(".pk")
    pickle.dump(lower_model, open(filename, "wb"))

    # MID Model
    # =========
    scaler = StandardScaler().fit(X_train)
    scaler.transform(X_train)
    scaler.transform(X_test)

    # mid_model.fit(X_train, y_train)
    objective = "regression"
    alpha = 0.9
    print("Prima di lgbm for MID model")
    mid_model = get_lgbm(
        X_train,
        X_test,
        y_train,
        y_test,
        cv,
        n_jobs,
        scoring,
        n_iter,
        objective,
        alpha,
        random_state,
    )
    print("Dopo lgbm for MID model")
    name = str(prop) + str("_mid")
    filename = str("pickle2-") + str(name) + str(".pk")
    pickle.dump(mid_model, open(filename, "wb"))

    # UPPER Model
    # ===========
    scaler = StandardScaler().fit(X_train)
    scaler.transform(X_train)
    scaler.transform(X_test)

    # upper_model.fit(X_train, y_train)
    objective = "quantile"
    alpha = UPPER_ALPHA
    print("Prima di lgbm for UPPER model")
    upper_model = get_lgbm(
        X_train,
        X_test,
        y_train,
        y_test,
        cv,
        n_jobs,
        scoring,
        n_iter,
        objective,
        alpha,
        random_state,
    )
    print("Dopo lgbm for UPPER model")
    name = str(prop) + str("_upper")
    filename = str("pickle2-") + str(name) + str(".pk")
    pickle.dump(upper_model, open(filename, "wb"))

    # PREDICTIONS and UQ
    lower = lower_model.predict(X_test)
    mid = mid_model.predict(X_test)
    upper = upper_model.predict(X_test)
    actual = y_test

    print("Model       mae     rmse")
    reg_sc = regr_scores(y_test, lower)
    info = OrderedDict()
    info["MAE_Lower"] = reg_sc["mae"]
    print("Lower:", round(reg_sc["mae"], 3), round(reg_sc["rmse"], 3))
    reg_sc = regr_scores(y_test, mid)
    info["MAE_Mid"] = reg_sc["mae"]
    print("Mid:", round(reg_sc["mae"], 3), round(reg_sc["rmse"], 3))
    reg_sc = regr_scores(y_test, upper)
    info["MAE_Upper"] = reg_sc["mae"]
    print("Upper:", round(reg_sc["mae"], 3), round(reg_sc["rmse"], 3))

    # Calculate the absolute error associated with prediction intervals
    # in_bounds = actual.between(left=lower, right=upper)
    if write_files:
        fout1 = open("Intervals.dat", "w")
        fout2 = open("Intervals1.dat", "w")
        line0 = "#    Jid      Observed       pred_Lower"
        line1 = "       pred_Mid        pred_Upper\n"
        line = line0 + line1
        fout1.write(line)
        line0 = "#    Jid      Observed       pred_Lower    AbsErr(Lower)"
        line1 = "     pred_Mid    AbsErr(Mid)     pred_Upper"
        line2 = "     AbsErr(Upper)    AbsErrInterval    Pred_inBounds\n"
        line = line0 + line1 + line2
        fout2.write(line)
        sum = 0.0
        count = 0
        MAE_err = 0.0
        for ii in range(len(actual)):
            true = float(actual[ii])
            llow = float(lower[ii])
            mmid = float(mid[ii])
            uupper = float(upper[ii])
            err = abs((uupper - llow) * 0.5)
            diff = true - mmid
            real_err = abs(diff)
            err_err = abs(real_err - err)
            MAE_err = MAE_err + err_err
            if abs(diff) < err:
                count = count + 1
                in_bounds = "True"
            else:
                in_bounds = "False"
            absolute_error_lower = abs(lower[ii] - actual[ii])
            absolute_error_mid = abs(mid[ii] - actual[ii])
            absolute_error_upper = abs(upper[ii] - actual[ii])
            absolute_error_interval = (
                absolute_error_lower + absolute_error_upper
            ) / 2.0

            line = (
                str(ii)
                + " "
                + jid[ii]
                + " "
                + str(actual[ii])
                + " "
                + str(lower[ii])
                + " "
                + str(mid[ii])
                + " "
                + str(upper[ii])
                + "\n"
            )
            sum = sum + float(absolute_error_interval)
            line2 = (
                str(ii)
                + " "
                + jid[ii]
                + " "
                + str(actual[ii])
                + " "
                + str(lower[ii])
                + " "
                + str(absolute_error_lower)
                + " "
                + str(mid[ii])
                + " "
                + str(absolute_error_mid)
                + " "
                + str(upper[ii])
                + " "
                + str(absolute_error_upper)
                + " "
                + str(absolute_error_interval)
                + " "
                + str(in_bounds)
                + "\n"
            )
            fout1.write(line)
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
    print("MAE predicted error (err=0.5*(High-Low))= " + str(MAE_error))
    info["MAE_Error"] = MAE_error
    return info
