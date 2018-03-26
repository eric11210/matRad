mask = zeros(ct.cubeDim);

tempMotionVecX = round(ct.motionVecX{2});
tempMotionVecY = round(ct.motionVecY{2});
tempMotionVecZ = round(ct.motionVecZ{2});

tempMotionVecX(tempMotionVecX < 1) = 1;
tempMotionVecY(tempMotionVecY < 1) = 1;
tempMotionVecZ(tempMotionVecZ < 1) = 1;

tempMotionVecX(tempMotionVecX > ct.cubeDim(1)) = ct.cubeDim(1);
tempMotionVecY(tempMotionVecY > ct.cubeDim(2)) = ct.cubeDim(2);
tempMotionVecZ(tempMotionVecZ > ct.cubeDim(3)) = ct.cubeDim(3);

V = sub2ind(ct.cubeDim,tempMotionVecY,tempMotionVecX,tempMotionVecZ);

mask(V) = 1;