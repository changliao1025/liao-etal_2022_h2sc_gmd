dArea = grc_pp%area(g) * 1000.0 * 1000.0 !unit is km, convert to m
dArea_half = dArea / 2.0
!this%rlen(g) is half the cell
dElevation_elm = ldomain%topo(g)
dChannel_length = dArea_half * this%gxr(g) !unit m
if(this%rlen(g) > dChannel_length) then
   dChannel_length = this%rlen(g)
end if
this%hlen(g) = dArea / dChannel_length
write(iulog,*) 'H2SC hillslope issue 1', this%hlen(g), dChannel_length

hlen_max = max(1250.0_r8, sqrt(dArea)/2.0 )
if (this%hlen(g) > hlen_max) then
   this%hlen(g) = hlen_max
end if