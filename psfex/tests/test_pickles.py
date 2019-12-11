import os
import io
import pickle
import numpy as np

from ..psfex_lib import PSFEx


def test_pickles():
    fname = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        'test_psfcat.psf'
    )

    psf = PSFEx(fname)
    im = psf.get_rec(10, 11)

    val = io.BytesIO()
    pickle.dump(psf, val)
    val.flush()
    val.seek(0)
    new_psf = pickle.load(val)
    new_im = new_psf.get_rec(10, 11)

    assert np.array_equal(new_im, im)

    # also make sure old object still works
    assert np.array_equal(im, psf.get_rec(10, 11))
