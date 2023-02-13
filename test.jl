using FITSIO


function main()
    f = FITS("path/to/file/filename.fits")

    fluxes = read(f[2], "PDCSAP_FLUX")
    dates = read(f[2], "TIME")
    println(dates)
end

main()