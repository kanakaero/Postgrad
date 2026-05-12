from dataclasses import dataclass, field
from typing import Tuple, List
import numpy as np


@dataclass
class Config:
    # Problem setting
    a: float = 0.1
    w: float = 100.0
    L: float = 1.0
    
    domain_spatio: np.ndarray = field(init=False)

    def __post_init__(self):
        self.domain_spatio = np.array([[0.], [self.L]], dtype=float)

    # Sampling strategy: "uniform" or "random"
    sampling: str = "random"
    
    # Uniform mesh sampling
    #dom_dis: List[int] = field(default_factory=lambda: [512])
    test_dis: List[int] = field(default_factory=lambda: [10000])
    # Random sampling
    n_dom: int = 128 # resample per epoch
    # n_bc: int = 200
    # n_ic: int = 200
    
    # Network
    n_hidden: int  = 2    
    w_neurons: int = 100 
    in_dim: int = 1      # x
    out_dim: int = 1      # u
    fourier_features: bool = True
    sigma: float = 1.0

    # Training / logging
    trials: int = 5 
    epochs: int = 40_000 
    disp: int = 10
    disp2: int = 1000     
    print_to_console: bool = True

    # Optimization
    use_decay_scheduler: bool = True # for Adam: (factor 0.9 every 1000 iterations)
    switch_to_lbfgs: bool = False 
    lbfgs_start_epoch: int = 0

    # ALM / CAPU params
    num_lambda: int = 1
    eta_vec: Tuple[float, ...] = (0.01,)

    # Paths
    out_dir: str = "logs"
    pic_dir: str = "pic"

    def make_method_name(self) -> str:
        return (
            f"psn_a{self.a}_w{int(self.w)}_L{self.L}"
            f"_batch_{int(self.n_dom)}"
            f"_fourier_nn{self.n_hidden}_{self.w_neurons}"
        )
    
def generate_configs():
    cases = [
        {"a": 0.1, "w": 100},
        {"a": 10.0, "w": 100}
    ]

    L_values = [1, 10, 100, 1000]

    configs = []
    for case in cases:
        for L in L_values:
            cfg = Config(
                a=case["a"],
                w=case["w"],
                L=L
            )
            configs.append(cfg)

    return configs