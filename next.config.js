/** @type {import('next').NextConfig} */
const nextConfig = {
    basePath: "/gene-platform",
    output: "export",  // <=== enables static exports
    reactStrictMode: true,
};

module.exports = nextConfig;