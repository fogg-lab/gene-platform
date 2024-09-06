import { LicenseInfo } from '@mui/x-license';

if (process.env.REACT_APP__MUI_KEY) {
  LicenseInfo.setLicenseKey(process.env.REACT_APP__MUI_KEY);
  console.log("License Info Set");
} else {
  console.warn("MUI license key not found in environment variables");
}
