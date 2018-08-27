#ifndef F11_H_
#define F11_H_

#include <genomes/GA1DArrayGenome.h>

long double f11_opt[1000] = {1.887532e+01,-1.906411e+01,-3.340450e+01,-4.840832e+01,-3.821750e+01,-2.716887e+01,3.100193e+01,-4.341158e+01,7.589020e+00,-3.750800e+01,-4.153965e+01,1.506957e+01,1.644970e+00,-1.863541e+01,-3.690158e+01,-2.296481e+01,2.922517e+00,1.515948e+01,8.648647e+00,3.520360e+01,5.350163e+00,-3.366402e+01,7.965197e+00,3.837632e+01,-2.898659e+01,-4.518140e+01,-4.741431e-01,6.828963e+00,-3.425629e+01,-4.498183e+01,5.345681e+00,-4.269115e+01,3.773348e+01,2.514779e+01,-2.299478e+01,1.805672e+01,-3.397367e+01,4.792290e+01,-4.814614e+01,1.479381e+01,-4.756185e+01,-3.397199e+01,2.448143e+01,-2.229243e+01,1.967725e+01,2.956010e+01,-3.451350e+00,-3.239824e+01,4.536837e+01,3.628956e+01,-2.054543e+01,2.528700e+01,-3.974154e+01,-4.383467e+01,2.345275e+01,1.786667e+01,-1.990498e+01,1.317384e+01,7.819404e+00,2.595462e+01,-2.817639e+01,-4.666188e+01,-3.460586e+01,-2.515605e+01,-3.304233e+01,-3.642485e+01,2.796758e+01,-1.923413e+01,4.994692e+01,-4.736564e+01,-3.161024e+00,-2.650909e+01,3.572907e+01,-8.003221e-01,-4.397094e+01,2.552173e+01,4.465887e+01,8.466440e+00,-6.095515e+00,-4.985939e+01,2.839725e+01,-2.690032e+00,2.811456e+00,3.661265e+01,6.076048e+00,-1.391555e+01,-1.877686e+01,3.198719e+01,-2.907958e+01,2.301187e+01,-1.834060e+01,3.215973e+01,1.321032e+00,-5.049543e-01,-2.666734e+01,1.158132e+01,-4.161836e+01,1.198768e+00,-3.612233e+01,4.967676e+01,-2.226582e+01,1.925234e+01,-3.017534e+01,4.040622e+01,3.595413e+01,-8.788558e+00,-4.673793e+01,2.502062e+01,6.404188e+00,4.642765e+00,-3.605763e+01,1.047043e+01,5.041910e+00,-3.232373e+01,-4.027905e+01,-1.782256e+01,4.058156e+01,9.394839e+00,-4.075116e+01,3.733672e+01,-1.761822e+01,3.432366e+01,1.334862e+01,-1.374252e+00,-3.403781e+01,-4.682630e+01,1.019691e+01,-1.682315e+01,5.546234e+00,-9.291341e+00,-5.901124e+00,-7.439305e+00,1.799342e+01,-4.562312e+01,3.381611e+01,3.418907e+01,1.719134e+01,6.985050e-01,-2.691117e+01,-4.939736e+01,9.041070e+00,-2.647982e+01,-4.920395e+01,-4.146276e+01,-2.952978e-01,2.802304e+01,2.958909e+01,-4.838773e+01,8.923357e-01,-4.525647e+01,2.528994e+01,-2.520997e+01,1.840923e+01,-4.974133e+01,-2.088827e+01,-2.210784e+01,-3.670999e+01,-3.437268e+01,3.126788e+01,-2.295669e+01,2.548160e+00,-1.481636e+01,-3.300872e+01,1.813949e+01,-3.869717e+01,-7.026575e+00,-1.930486e+01,-1.017086e+01,-1.245818e+01,2.193656e+01,-1.741543e+01,2.500956e+01,1.137915e-01,3.707020e+01,-3.224460e+01,-2.084010e+01,1.279794e+01,2.045145e+01,-4.977200e+01,1.079185e+01,3.115626e+01,-7.261440e+00,3.945303e+01,-2.454991e+01,-6.083330e+00,1.794258e+01,3.073471e+01,-2.462946e+01,-4.149694e+01,3.774959e+01,-1.529533e+01,4.728833e+01,1.275971e+01,-2.252183e+01,-4.799545e+01,-2.442176e+01,3.279052e+01,4.842036e+01,3.500368e+00,3.748993e+01,1.850586e+01,2.784867e+01,9.026225e+00,3.299667e+01,1.641658e+01,1.816176e+01,-9.177620e+00,-1.951185e+01,-3.602836e+01,1.691951e+01,2.621134e+01,-4.639620e+01,-1.693995e+01,3.058611e+01,3.046602e+00,-4.253892e+01,1.640818e+01,4.971472e+01,-2.625811e+01,1.976268e+01,-4.808116e+01,2.745219e+01,3.450348e+01,-4.549077e+01,7.192325e-01,1.100263e+01,-5.661636e+00,-3.679276e+01,3.294121e+01,2.024530e+01,6.311754e+00,2.657141e+01,-1.580554e+01,-2.441378e+01,-2.917622e+01,-2.433787e+01,-1.975400e+01,2.998845e+00,-1.858891e+01,3.692511e+01,-3.702805e+01,3.703302e+00,-2.101376e+01,-3.315787e+01,4.608837e+01,4.266524e+01,-3.164042e+01,-4.985673e+01,-7.012149e+00,2.729015e+01,-3.936018e+01,1.242190e+01,-3.514933e+00,-2.352390e+01,3.530878e+01,-6.818319e+00,1.916810e+00,9.349743e+00,2.891005e+00,4.137915e+01,-1.552180e+01,5.131123e+00,-4.823620e+01,-4.606989e+01,1.004916e+01,-7.463954e+00,-4.727222e+01,-4.397388e+01,-2.452148e+01,4.246483e+01,1.032807e+00,3.802185e+01,2.492280e+00,-2.726305e+00,-2.320402e+01,-2.555842e+01,3.053696e+01,2.423886e+01,-4.816449e+01,-2.371661e+01,4.697462e+01,2.190729e+01,3.517454e+00,-2.768089e+01,-2.485074e+01,-2.569007e+01,3.845727e+01,1.899002e+01,-2.418130e+01,4.671832e+01,-3.387269e+01,-3.117013e+01,4.432485e+01,3.431091e+01,-4.160253e+01,-3.701950e+01,1.380134e+00,3.256118e+00,3.511257e+01,-4.166794e+01,4.362809e+01,-4.967984e+01,1.005182e+01,-6.632051e+00,-1.368292e+01,1.314863e+01,-2.659956e+01,4.458856e+01,4.593362e+01,8.654529e+00,3.983999e+01,-4.275081e+01,-2.132565e+01,-4.688008e+01,-2.575239e+01,3.935569e+01,4.947831e+01,5.900004e+00,-3.902671e+01,2.012107e+01,3.571731e+01,4.482406e+00,-4.360513e+01,1.229474e+01,-4.348048e+01,-1.787087e+01,-2.402584e+01,-2.398424e+01,4.170561e+01,-4.096838e+01,-4.456630e+01,3.565386e+01,1.448472e+01,-9.937257e+00,-4.360807e+01,-2.477357e+01,4.104555e+01,-2.435048e+01,1.530892e+01,-1.598298e+01,-1.718077e+00,1.514268e+01,-1.899142e+01,-4.944295e+00,3.303995e+01,2.984230e+01,9.549386e-01,4.192339e+01,-2.149918e+01,-4.475929e+01,1.729484e+01,-1.974238e+01,-2.191051e+01,2.387332e+01,3.380911e+01,1.645131e+01,4.028900e+01,1.131046e+01,8.976506e+00,4.699940e+01,-1.163972e+01,3.646882e+01,-4.112811e+00,1.222289e+01,-4.048017e+01,3.560169e+00,3.848570e+01,3.101131e+01,3.758762e+00,-4.274451e+01,2.121823e+01,-2.872329e+01,-1.283211e+01,2.512727e+00,1.150135e+01,-2.722069e+01,1.992696e+01,-1.538595e+01,8.320367e+00,2.433395e+01,5.600994e+00,3.675340e+01,1.722566e+01,9.396660e+00,-1.232429e+01,-3.918665e+01,-2.948195e+01,-7.731872e+00,1.175092e+01,-3.284402e+01,-8.568818e+00,2.712069e+01,-1.050586e+01,1.298995e+01,2.053352e+01,-1.212773e+00,3.613480e+01,7.822345e+00,4.259046e+01,4.894051e+01,-3.816743e+00,-4.354434e+01,1.540569e+01,-1.791625e+01,3.043234e+01,2.110605e+01,8.006092e+00,-4.497427e+01,4.206947e+01,3.628535e+01,-4.649788e+01,2.994594e+01,-1.733840e+01,-3.534099e+01,4.930255e+01,-5.392178e+00,4.557473e+00,2.001541e+00,4.216400e+01,-2.836014e+01,-3.192556e+01,-1.805350e+01,-2.464550e+00,-1.990414e+01,6.969434e+00,4.354084e+01,-6.429117e+00,-4.377907e+01,3.185302e+01,1.082924e+01,3.133595e+01,-4.659431e+00,-3.594699e+01,-3.244613e+01,1.554911e+01,-1.764749e+01,4.738385e+01,3.615427e+01,3.117895e+01,4.962914e+01,7.754700e+00,4.502412e+01,4.994734e+01,-5.406603e+00,-3.545625e+01,-2.535415e+00,-1.575960e+01,-1.495501e+01,2.463058e+01,4.969441e+01,-5.353454e-01,7.454151e+00,1.329582e+01,-3.833318e+01,-2.052372e+01,6.530234e+00,5.416827e+00,-4.630587e+01,-3.842169e+01,-1.223255e+01,2.385778e+01,-5.422569e+00,3.985036e+01,-3.279528e+01,-2.396338e+01,4.099891e+01,-1.688617e+01,2.884990e+01,2.123378e+01,2.430230e+01,1.651714e+01,-3.792381e+01,3.229950e+01,2.113813e+01,3.292483e+01,-3.359056e+00,-1.802647e+01,1.268646e+01,-1.044130e+01,-9.670180e+00,3.948160e+01,4.459991e+01,-8.635902e+00,3.715969e+01,1.731151e+01,1.865306e+01,3.258044e+01,-4.966766e+01,-1.104114e+01,-4.214005e+01,5.236581e+00,5.456182e+00,-6.791849e+00,1.874234e+00,3.479255e+01,4.259466e+01,-2.245510e+00,4.740373e+01,8.372956e-01,-1.312020e+01,-3.107657e+01,7.814642e+00,-1.428809e+01,4.288064e+01,3.322636e+01,3.983103e+01,4.163475e+01,-6.739330e+00,4.624173e+01,-2.551010e+01,-2.271342e+01,-3.026624e+01,1.852666e+00,3.733294e+01,-1.445307e+01,-4.644200e+01,2.427611e+01,4.644662e+01,2.172564e+01,2.333763e+01,-1.777116e+01,-2.806785e+01,3.489745e+01,-3.207654e+01,3.239754e+01,3.525052e+01,4.283779e+01,3.135836e+01,-5.546514e+00,-4.267568e+00,1.197710e+01,3.044256e+01,1.285914e+01,-4.084038e+01,8.379328e+00,1.821554e+01,-4.636259e+01,7.447848e+00,3.106523e+01,-1.682245e+01,-3.479241e+01,6.717622e+00,-1.896887e+01,1.121018e+01,4.581149e+01,-4.639214e+01,2.475102e+01,4.489318e+01,1.292539e+01,6.234726e+00,2.850243e+01,4.470117e+01,1.273184e+01,-8.334792e+00,-3.717832e+01,-3.289108e+01,4.419614e+01,-4.693274e+01,2.288764e+01,2.559721e+01,-1.856595e+01,-1.883337e+00,2.925927e+01,-2.365891e+01,4.306159e+01,2.257813e+01,3.814313e+01,-1.555639e+01,1.650467e+01,1.721922e+01,2.824488e+01,3.966857e+01,3.324744e+00,2.929232e+01,2.531179e+01,-2.534918e+01,4.835622e+01,3.692595e+01,-3.356024e+01,-6.781345e+00,-1.500298e+00,-8.155247e+00,9.041210e+00,-4.879010e+01,2.209047e+01,2.881433e+01,-1.178369e+01,-1.602934e+01,-4.589888e+01,2.721186e+01,3.127320e+01,-4.034922e+01,8.314905e+00,-1.783194e+01,4.460558e+00,-3.245706e+01,-1.245244e+01,-7.751759e+00,4.333651e+01,2.186835e+01,-1.384048e+01,3.170904e+01,-3.889423e+01,-3.914296e+01,1.142712e+01,-1.506740e+00,7.666328e+00,4.138139e+01,4.672574e+01,-4.789398e+00,3.672231e+01,-3.394594e+01,9.777669e-01,4.740891e+01,3.929449e+01,2.445902e+01,4.253612e+01,-4.541865e+01,-2.591471e+01,-1.133091e+01,7.668289e+00,-4.553972e+00,3.381261e+01,-3.609741e+01,3.789258e+01,2.020917e+01,6.349848e+00,3.733910e+01,3.571437e+01,1.487798e+01,-4.696789e+01,1.805784e+01,1.270397e+01,-3.899128e+01,-1.942068e+01,4.351549e+01,2.330262e+01,3.150317e+01,4.153755e+01,-4.959385e+01,-1.968313e+01,-1.540205e+01,-2.314660e+01,2.572872e+01,-2.538980e+01,1.872715e+01,6.818179e+00,-3.298211e+01,3.780533e+01,-3.400644e+01,-4.304422e+01,-4.095900e+01,-1.419397e+01,2.784503e+01,-1.425090e+00,6.020728e+00,1.253647e+01,-1.523581e+01,-4.571990e+01,2.631596e+01,1.457869e+01,-3.356955e+00,2.511040e+01,-1.633983e+01,4.601485e+01,2.593067e+01,2.013116e+01,1.833192e+01,-7.499807e+00,6.761738e+00,2.639522e+01,-3.083919e+01,-1.431399e+01,3.380729e+01,-2.551066e+01,-2.970617e+01,1.656504e+01,-4.712139e+01,-1.446875e+01,-1.268898e+01,-1.423711e+01,-4.338342e+01,1.022296e+01,-1.106509e+01,-4.758720e+01,-4.539862e+01,-2.437919e+01,4.619236e+00,9.179020e+00,-3.119814e+01,3.412605e+01,3.798347e+01,-3.524575e+01,-2.674213e+01,6.203494e+00,1.424523e+01,4.480368e+01,1.321151e+01,-4.204524e+01,-2.494724e+01,6.547880e+00,-4.396688e+01,-1.790757e+01,1.386387e+01,-2.270964e+01,1.775183e+01,3.314380e+00,-1.784987e+01,3.297153e+00,6.071216e-02,-2.156122e+01,-8.621197e+00,2.324463e+01,-1.319233e+01,-7.195897e+00,-1.387676e+01,-4.035272e+01,3.403921e+01,-2.286972e+01,-1.903568e+01,4.414894e+01,-2.195742e+01,-3.805056e+01,-7.420118e+00,4.931837e+01,3.354371e+01,6.215119e+00,4.989944e+01,1.375799e+01,-2.495620e+01,-4.950310e+01,-1.452043e+01,-2.440538e+01,-2.975316e+00,-2.703638e+01,1.353594e-01,1.844130e+01,2.643850e+01,2.339897e+01,-3.016946e+01,2.883225e+01,-3.496593e+01,-4.880746e+01,-4.552341e+01,-1.915416e+01,3.622485e+01,3.127783e+01,-1.437562e+01,2.895676e+01,6.565526e+00,4.477973e+01,2.748146e+01,-4.811015e+01,2.109205e+01,4.293946e+01,3.192948e+01,-4.958405e+01,6.850530e+00,1.548482e+01,2.616512e+01,-4.248038e+01,4.976527e+01,4.024964e+01,4.983572e+01,3.388236e+01,-3.545121e+01,3.964126e+01,3.104870e+01,4.682609e-01,2.539260e+01,3.858408e-02,-4.619376e+00,1.757200e+01,-3.160401e+01,9.673961e+00,-6.390322e+00,-4.166332e+01,3.569742e+01,-4.785554e+01,-3.198410e+01,1.228409e+01,3.011568e+01,-1.290172e+01,4.002724e+01,-4.556808e+01,-1.712944e+01,-4.593502e+01,2.741102e+01,-1.986086e+01,-2.768594e+01,-3.574546e+01,2.356871e+01,3.465166e+01,4.161892e+01,3.932236e+01,-3.464886e+01,-4.751521e+01,3.915570e+01,-2.879234e+01,3.766045e+00,1.553622e+01,-4.639172e+01,-2.388523e+01,4.703568e+01,1.497798e+01,1.200161e+01,-4.555310e+01,-1.731025e+01,2.533308e+01,-3.893554e+01,-4.098575e+01,4.635881e+01,7.457092e+00,-8.944015e+00,1.803235e+01,-2.500648e+01,2.606526e+01,-6.260075e+00,1.601408e+01,-4.400917e+01,3.377172e+01,1.011239e+00,-2.727390e+01,-2.690025e+01,1.553482e+01,3.865138e+01,3.109814e+01,4.079374e+01,4.348468e+01,2.525192e+00,-1.028325e+00,-4.733301e+01,-4.105599e-01,-1.481748e+01,1.077658e+01,-2.376142e+01,1.496691e+01,2.865943e+01,-8.054830e+00,-5.930675e+00,-4.852386e+01,4.119849e+01,1.618886e+01,-1.246224e+01,9.373691e+00,-4.879584e+01,-4.530591e+01,7.907076e+00,-3.130107e+01,4.195168e+01,-4.908112e+01,3.427520e+01,4.857316e+01,-6.885123e+00,2.199314e+01,2.504884e+00,3.614096e+01,-2.732047e+00,-2.371745e+01,-1.341557e+01,-3.676167e+01,2.142509e+01,-4.320150e+01,3.338588e+01,6.646196e+00,-3.176562e+01,-4.217913e+01,2.383803e+01,1.211232e+00,2.206764e+01,-4.516249e+01,2.747712e+01,4.391058e+01,-2.071083e+01,-3.853780e+01,3.422464e+01,2.603739e+01,3.834579e+01,1.991450e+01,4.510535e+01,-2.453380e+01,-3.765211e+01,-4.795904e+01,-4.651805e+01,-3.466776e+01,4.079682e+01,-1.673632e+01,4.541459e+01,2.251749e+01,4.330093e+01,-1.488652e+01,3.328315e-01,-3.709429e+01,-4.064150e+01,4.834193e+01,3.821190e+01,3.369875e+01,-4.090970e+01,1.330969e+01,-4.785806e+01,-1.845405e+01,-2.462932e+01,3.002199e+01,4.560183e+01,2.462904e+01,3.943945e+01,1.030741e+01,2.577998e+01,-1.152082e+01,1.302847e+01,-1.009958e+01,-4.432877e+01,4.771409e+01,1.593523e+01,4.546809e+01,4.289857e+01,4.189566e+01,-2.537579e+01,-4.344239e+01,-2.083505e+01,-1.077532e+01,3.713238e+01,-4.055530e+00,3.045671e+01,-2.566248e+01,1.712860e+01,-4.693176e+01,-2.488911e+01,4.509485e+01,-4.882595e+01,4.636917e+01,-4.783397e+01,-5.299044e+00,-7.566192e+00,-3.141606e+01,3.301866e+01,1.044732e+01,4.783936e+00,-6.533595e+00,-3.141760e+01,-4.079332e+01,-3.036693e+01,5.296803e+00,2.825636e+01,-3.921172e+01,1.986982e+01,3.624733e+00,1.533007e+01,-4.040965e+00,-4.962746e+01,7.800077e+00,4.876930e+00,2.208291e+01,-4.337866e+01,-3.119758e+01,-3.117251e+01,3.974728e+01,-4.530633e+01,1.587823e+01,2.587760e+01,5.988936e+00,4.217576e+01,-3.942460e+01,-2.341788e+01,1.944281e+01,-1.768853e+01,-2.780190e+01,4.297055e+01,-2.511754e+01,3.263786e+01,4.745765e+01,3.113091e+01,-1.998368e+01,-4.701915e+01,3.499996e+01,2.567522e+01,-3.975344e+01,2.971065e+01,-2.930325e+01,-4.894275e+01,-4.636315e+01,1.142488e+01,4.248528e+01,-3.221827e+01,-1.677399e+01,-5.524526e+00,-1.891565e+01,-2.739827e+01,2.933700e+01,-1.946255e+01,-1.272750e+01,3.279507e+00,2.389993e+01,5.717937e+00,-3.206687e+01,-4.636371e+01,2.658023e+01,8.831273e+00,-4.906124e+01,-8.061973e+00};

// F11 also known as Schaffer

const long double f_Schaffer_BIAS = 0.0;

long double f_Schaffer(GAGenome& g) {

   GA1DArrayAlleleGenome<long double>& genome = DYN_CAST (GA1DArrayAlleleGenome<long double>&, g);

   unsigned dim = genome.length ();

   long double sum, aux, aux2;
   long double currentGen, nextGen;

   sum = 0.0;
   currentGen = genome.gene (0) - f11_opt[0];
   currentGen = currentGen * currentGen;

   for (unsigned i = 1; i < dim; i++) {
      nextGen = genome.gene (i) - f11_opt[i];
      nextGen = nextGen * nextGen;
      aux = currentGen + nextGen;
      currentGen = nextGen;
      aux2 = sin (50.0 * pow (aux, 0.1));
      sum += pow (aux, 0.25) * (aux2 * aux2 + 1.0);
   }

   return sum;

}

#endif /* F11_H_ */