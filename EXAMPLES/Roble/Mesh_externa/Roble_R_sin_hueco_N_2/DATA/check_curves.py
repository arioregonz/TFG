
import gmsh

gmsh.initialize()

gmsh.open("three.geo")

gmsh.model.occ.synchronize()



print("\n=== CURVAS Y SUS POSICIONES ===")

for c in gmsh.model.getEntities(1):

    bbox = gmsh.model.getBoundingBox(1, c[1])

    xmin,ymin,zmin,xmax,ymax,zmax = [round(x,1) for x in bbox]

    

    if ymin==0 and ymax==0:

        label = "<<< BOTTOM"

    elif xmax==100 and xmin==100:

        label = "<<< RIGHT"

    elif ymin==100 and ymax==100:

        label = "<<< TOP"

    elif xmin==0 and xmax==0:

        label = "<<< LEFT"

    else:

        label = "(interior/curva)"

    

    print(f"Curve {c[1]:3d}: x=[{xmin},{xmax}] y=[{ymin},{ymax}] {label}")



gmsh.finalize()

