# Change soil layer thickness



```
  dzmm(c,j) = dz(c,j)*1.e3_r8
  dz                 =>    col_pp%dz  

```

```

  subroutine initVertical(bounds, snow_depth, thick_wall, thick_roof)

```

```
    dzsoi(1) = 0.5_r8*(zsoi(1)+zsoi(2)) !thickness b/n two interfaces
    do j = 2,nlevgrnd-1
       dzsoi(j)= 0.5_r8*(zsoi(j+1)-zsoi(j-1))
    enddo
    dzsoi(nlevgrnd) = zsoi(nlevgrnd)-zsoi(nlevgrnd-1)
```

```
    real(r8), allocatable :: dzsoi(:)        !soil dz (thickness)
```

```
    dzsoi(nlevgrnd) = zsoi(nlevgrnd)-zsoi(nlevgrnd-1)
    !h2sc change last layer thickness
    dzsoi(nlevgrnd) = 100.0
```