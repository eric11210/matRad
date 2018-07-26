function dicomInfo = matRad_dicomHeader(xcatLog)

StudyInstanceUID = dicomuid;
FrameOfReferenceUID = dicomuid;

dicomInfo.PixelSpacing = [xcatLog.resolution.x; xcatLog.resolution.y];
dicomInfo.SliceThickness = xcatLog.resolution.z;
dicomInfo.ImageOrientationPatient = [1; 0; 0; 0; 1; 0];
dicomInfo.PatientPosition = 'HFS';
dicomInfo.SpacingBetweenSlices = -xcatLog.resolution.z;
dicomInfo.DataCollectionDiameter = 600;
dicomInfo.ReconstructionDiameter = 435.2;

dicomInfo.StudyInstanceUID = StudyInstanceUID;
dicomInfo.SOPInstanceUID = dicomuid;
dicomInfo.SOPClassUID = '1.2.840.10008.5.1.4.1.1.2';

dicomInfo.MediaStorageSOPInstanceUID = dicomInfo.SOPInstanceUID;
dicomInfo.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.2';

dicomInfo.TransferSyntaxUID = '1.2.840.10008.1.2.1';
dicomInfo.ImplementationClassUID = '1.2.276.0.7230010.3.0.3.6.0';

dicomInfo.StudyDate = '20160603';
dicomInfo.AcquisitionDate = '20160603';


dicomInfo.PatientID = '4815162342';
dicomInfo.PatientBirthDate = '18670701';
dicomInfo.PatientSex = 'M';
dicomInfo.PatientAge = '148Y';

dicomInfo.StudyID = '27570';

dicomInfo.FrameOfReferenceUID = FrameOfReferenceUID;

dicomInfo.Modality = 'CT';

dicomInfo.PhotometricInterpretation = 'MONOCHROME2';
dicomInfo.BitsAllocated = 16;
dicomInfo.BitsStored = 16;
dicomInfo.SmallestImagePixelValue = 0;
dicomInfo.LargestImagePixelValue = 65535;


dicomInfo.Manufacturer = 'XCAT';
dicomInfo.ManufacturerModelName = 'XCAT';
dicomInfo.ConvolutionKernel = 'XCAT';

end

