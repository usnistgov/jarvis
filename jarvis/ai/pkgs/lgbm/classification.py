"""Module for LightGBM classification."""

from jarvis.ai.pkgs.sklearn.classification import classification as sk_cl


def classification(
    X=[],
    Y=[],
    tol=100,
    plot=False,
    preprocess=False,
    models=[],
    model_name="my_model",
    save_model=False,
):
    """Provide function for classification models."""
    info = sk_cl(
        X=X,
        Y=Y,
        tol=tol,
        preprocess=preprocess,
        models=models,
        model_name=model_name,
        save_model=save_model,
    )
    return info


"""
if __name__ == "__main__":
    from jarvis.ai.pkgs.utils import get_ml_data, binary_class_dat
    property = "exfoliation_energy"
    tol=100
    #property = 'optb88vdw_bandgap'
    #tol=0.05
    X, Y, jid = get_ml_data(dataset="cfid_3d", ml_property=property)
    #X_class, Y_class = binary_class_dat(X=X, Y=Y, tol=tol)
    print ('lennnn',len(X))
    #models = [LGBMClassifier(n_estimators=1000,max_depth=50,num_leaves=100)]
    models = [LGBMClassifier()]
    info = classification(X=X, Y=Y, models=models,
    preprocess=False, save_model=False, tol=tol)
    print (info['LGBMClassifier']['roc_auc'][0])
"""
