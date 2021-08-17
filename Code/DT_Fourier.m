function X = DT_Fourier(x,N0,w)
    X = 0;
    for x_index = 1 : size(x,2) %iterating through every element in x(t)
        X = X + x(x_index) .* exp((-1.*1j) .* w .* (x_index-N0));  %calculating Discrete Time Fourier transform for every w stored in w
    end 
end

