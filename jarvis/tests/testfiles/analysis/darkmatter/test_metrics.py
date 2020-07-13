from jarvis.analysis.darkmatter.metrics import peskin_schroeder_constant
import numpy as np

def test_peskin_schroeder():
    epsilon = [187.5, 9.8, 90.9]
    fermi_velocities = [0.0029, 0.0050, 0.0021]
    val = peskin_schroeder_constant(
        bandgap=0.035, fermi_velocities=fermi_velocities, epsilon=epsilon
    )
    assert round(val[0],5)==round(1.83908046,5)
    assert round(val[1],5)==round(20.40816327,5)
    assert round(val[2],5)==round(5.2386191,5)
