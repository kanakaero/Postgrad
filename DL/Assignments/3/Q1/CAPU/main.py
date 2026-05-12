import torch

from config import Config, generate_configs
from core import Trainer
from utils import _torch_device, _ensure_dirs

def main():
    torch.set_default_dtype(torch.float64)
    configs = generate_configs()

    for cfg in configs:
        print("\n==============================")
        print(f"Running: a={cfg.a}, w={cfg.w}, L={cfg.L}")
        print("==============================")

        trainer = Trainer(cfg)
        trainer.run()

if __name__ == "__main__":
    main()
