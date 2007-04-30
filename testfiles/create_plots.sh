#!/bin/sh
cd convergence

echo "Linear Interpolant C2F"
cd linear_interpolant_C2F
gnuplot interpolation
cd ..

echo "Linear Interpolant F2C"
cd linear_interpolant_F2C
gnuplot interpolation
cd ..

echo "Gradient C2F"
cd gradient_C2F
gnuplot gradient
cd ..

echo "Gradient F2C"
cd gradient_F2C
gnuplot gradient
cd ..

echo "Divergence C2F"
cd divergence_C2F
gnuplot divergence
cd ..

echo "Divergence F2C"
cd divergence_F2C
gnuplot divergence
cd ..

echo "Scratch Laplacian"
cd scratch_laplacian
gnuplot laplacian
cd ..

echo "Mixed Derivatives"
cd mixed_derivatives
gnuplot mixed_derivatives
cd ..

cd ..
pdflatex report.tex


