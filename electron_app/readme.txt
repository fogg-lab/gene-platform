Note: This is a rough readme and stand-in for the future readme.md

Install and run steps:
1. install dependencies: npm install
2. start webpack dev server: npm start
3. run electron application: npm run electron

This should launch the electron/react window.

Project Structure

electron_app/
│
├── src/
│   ├── index.js
│   ├── App.js //handles React (all UI goes here)
│   └── styles.css (not implemented)
│
├── public/
│   └── index.html //root React element
│
├── dist/
│   └── (output files after build)
│
├── main.js //handles electron
├── package.json
├── webpack.config.js
├── .babelrc
└── README.md

