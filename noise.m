function nu = noise(t)
    std_tacho = 0.01; % tachometer standard deviation
    std_GPS = 0.05; % gps standard deviation

    nu = 0*[
        std_tacho*randn(1);
        std_GPS*randn(1)
        ];
end

