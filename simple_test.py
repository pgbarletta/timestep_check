import concurrent.futures as cf
from typing import Optional, Any

import MDAnalysis as mda  # type: ignore
import numpy as np
from MDAnalysis.transformations.base import TransformationBase  # type: ignore
from numpy.typing import NDArray

from MDAnalysis.tests.datafiles import TPR, XTC


class one_transform(TransformationBase):
    def __init__(
        self,
        ag: mda.AtomGroup,
        max_threads: Optional[int] = None,
        parallelizable: Optional[bool] = True,
    ):
        super().__init__(max_threads=max_threads, parallelizable=parallelizable)
        self.ag = ag

    def _transform(self, ts: Any):  # type: ignore
        print(f"Frame inside the transform: {ts.frame}")
        return ts


u = mda.Universe(TPR, XTC)
f = one_transform(u.atoms)
# Serial
# for idx, ts in enumerate(u.trajectory):
#     print(f"Frame input: {ts.frame}")
#     analysis_func(ts)

# Parallel
with cf.ProcessPoolExecutor(max_workers=4) as ex:
    futures = []
    for idx, ts in enumerate(u.trajectory):
        print(f"Frame input: {ts.frame}")
        futures.append(ex.submit(f, ts))
    for futu in cf.as_completed(futures):
        if futu.exception():
            print(f"Exception:  {futu.exception()} ")
        ts = futu.result()

