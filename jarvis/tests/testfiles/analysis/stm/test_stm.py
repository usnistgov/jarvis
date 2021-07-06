from jarvis.analysis.stm.tersoff_hamann import TersoffHamannSTM
import matplotlib.pyplot as plt
import os
from jarvis.db.figshare import make_stm_from_prev_parchg
from jarvis.core.image import Image

name = os.path.join(os.path.dirname(__file__), "PARCHG")
from jarvis.core.image import Image
from io import BytesIO


def test_th_stm():
    plt.switch_backend("agg")
    TH_STM1 = TersoffHamannSTM(chg_name=name, min_size=5)
    # byte_io = BytesIO()
    # t1 = TH_STM1.constant_height(filename=byte_io)

    t1 = TH_STM1.constant_height(filename="test.png")
    im = Image.from_file("test.png")
    plt.imshow(
        im.fourier_transform2D(use_crop=True, zoom_factor=10)
        .rotate(angle=0)
        .black_and_white(threshold=0.05)
        .values,
        cmap="Greys",
    )
    # p = byte_io.getvalue()  # .decode('UTF-8')
    # print ('p',p)
    t1 = TH_STM1.constant_height()
    TH_STM2 = TersoffHamannSTM(chg_name=name, min_size=5)
    t2 = TH_STM2.constant_current()
    t2 = TH_STM2.constant_current(pc=5)
    from jarvis.db.figshare import make_stm_from_prev_parchg

    make_stm_from_prev_parchg()
    Image.get_blob_angles(filename="stm_image.png")
    im = Image.crop_from_center()
    ims = Image.augment_image()

    cmd = "rm *.png"
    os.system(cmd)


# test_th_stm()
