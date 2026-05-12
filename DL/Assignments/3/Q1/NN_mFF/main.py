# -*- coding: utf-8 -*-

import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import seaborn as sns
from models_tf import Sampler, NN_mFF

import os
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
os.makedirs("plots", exist_ok=True)

if __name__ == '__main__':

    model_list = [
        ("NN_mFF", NN_mFF)
    ]

    cases = [
        {"a": 0.1, "w": 100},
        {"a": 10,  "w": 100}
    ]

    L_values = [1, 10, 100, 1000]  

    # Exact solution + noise
    def u(x, a, w):
        noise = 0.05 * np.random.randn(*x.shape)
        return np.sin(2*np.pi*x) + a*np.sin(2*np.pi*w*x) + noise

    # Exact PDE residual
    def u_xx(x, a, w):
        return -4*np.pi**2*np.sin(2*np.pi*x) - 4*a*(w**2)*(np.pi**2)*np.sin(2*np.pi*w*x)

    for case in cases:
        for L in L_values:
            for model_name, ModelClass in model_list:

                print(f"\n===== {model_name} | a={case['a']} | w={case['w']} | L={L} =====")

                a = case['a']
                w = case['w']

                # Define computational domain
                bc1_coords = np.array([[0.0], [0.0]])
                bc2_coords = np.array([[L], [L]])
                dom_coords = np.array([[0.0], [L]])

                # Create boundary sampler
                bc1 = Sampler(1, bc1_coords, lambda x: u(x, a, w), name='Dirichlet BC1')
                bc2 = Sampler(1, bc2_coords, lambda x: u(x, a, w), name='Dirichlet BC2')

                bcs_samplers = [bc1, bc2]

                # Create residual sampler
                res_samplers = Sampler(1, dom_coords, lambda x: u_xx(x, a, w), name='Forcing')

                # Define model
                if model_name == "NN":
                    layers = [1, 100, 100, 1]
                else:
                    layers = [100, 100, 1]
    
                # Hyper-parameter for Fourier features
                sigma = 10
                
                # NN: Vanilla MLP
                # NN_FF : Vanilla Fourier feature network
                # NN_mFF : Multi-scale Fourier feature network
                model = ModelClass(layers, bcs_samplers, res_samplers, u, a, w, sigma)

                # Train model
                model.train(nIter=40000, batch_size=128, log_NTK=False, log_weights=False)

                # Create test data
                nn = 10000
                X_star = np.linspace(dom_coords[0, 0], dom_coords[1, 0], nn)[:, None]
                u_star = u(X_star, a, w)
                r_star = u_xx(X_star, a, w)

                # Predictions
                u_pred = model.predict_u(X_star)
                r_pred = model.predict_r(X_star)
                error_u = np.linalg.norm(u_star - u_pred, 2) / np.linalg.norm(u_star, 2)
                error_r = np.linalg.norm(r_star - r_pred, 2) / np.linalg.norm(r_star, 2)

                print('Relative L2 error_u: {:.2e}'.format(error_u))
                print('Relative L2 error_r: {:.2e}'.format(error_r))
                        
                loss_bcs = model.loss_bcs_log
                loss_res = model.loss_res_log
                l2_error = model.l2_error_log
    
                # Plot
                fig = plt.figure(figsize=(18, 5))
                with sns.axes_style("darkgrid"):
                    plt.subplot(1, 3, 1)
                    plt.plot(X_star, u_star, label='Exact')
                    plt.plot(X_star, u_pred, '--', label='Predicted')
                    plt.xlabel('$x$')
                    plt.ylabel('$y$')
                    plt.tight_layout()

                    plt.subplot(1, 3, 2)
                    plt.plot(X_star, u_star - u_pred, label='Error')
                    plt.xlabel('$x$')
                    plt.ylabel('Point-wise error')
                    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
                    plt.tight_layout()

                    plt.subplot(1, 3, 3)
                    iters = 100 * np.arange(len(loss_res))

                    plt.plot(iters, loss_res, label='$\mathcal{L}_{r}$', linewidth=2)
                    plt.plot(iters, loss_bcs, label='$\mathcal{L}_{b}$', linewidth=2)
                    plt.plot(iters, l2_error, label=r'$L^2$ error', linewidth=2)

                    plt.yscale('log')
                    plt.xlabel('iterations')
                    plt.legend(loc='upper right', bbox_to_anchor=(1.0, 0.9), fontsize=20)

                    fig.text(
                                0.5, -0.05,
                                f"Multiscale Fourier Feature Neural Network with a={a}, w={w}, L={L}",
                                ha='center', fontsize=14
                            )

                    plt.tight_layout()
                    plt.savefig(f"plots/{model_name}_a1{a}_a2{w}_L{L}.png", dpi=600, bbox_inches='tight')