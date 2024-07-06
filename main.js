const { app, BrowserWindow } = require('electron');
const path = require('path');
const { default: installExtension, REACT_DEVELOPER_TOOLS } = require('electron-extension-installer');

let mainWindow;

async function createWindow() {
  const isDev = (await import('electron-is-dev')).default;

  mainWindow = new BrowserWindow({
    width: 1200,
    height: 900,
    webPreferences: {
      nodeIntegration: true,
      contextIsolation: false,
    }
  });

  const url = isDev
    ? 'http://localhost:9000'
    : `file://${path.join(__dirname, 'dist/index.html')}`;

  console.log('Loading URL:', url);

  mainWindow.loadURL(url).catch(err => {
    console.error('Failed to load URL:', err);
  });

  mainWindow.on('closed', () => {
    mainWindow = null;
  });
}

app.on("ready", async () => {
  await installExtension(REACT_DEVELOPER_TOOLS, {
    loadExtensionOptions: {
      allowFileAccess: true,
    },
  });
  createWindow();
});

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

app.on('activate', () => {
  if (mainWindow === null) {
    createWindow();
  }
});
